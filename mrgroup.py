import sqlite3
import lmdb
import pandas as pd
import numpy as np

from operator import itemgetter
from scipy.signal import find_peaks

from database import get_rid2chrom_conn
from utils import return_sm_matrix, return_sm_vector, rbf_kernel_numba_matrix, rbf_kernel_numba_vector, return_reps_vector
from chromatographic import organize_chroms, return_xics_by_fid

class PT_BM():
    def __init__(self, fid2repss, fid2nrts):
        self.fid2repss = fid2repss
        self.fid2nrts = fid2nrts
        self.rid_list = []
        self.fid_list = []
        self.rid_u_list = []
        self.fid_u_list = []
        self.mrss_list = []
        self.mrs_list = []
        
        self.ptp_peaks = None
        self.ptp_intens = None

def return_bm_from_pt(p_table: pd.DataFrame, nrt_intervel: float, nrt_width: float):
    '''Extract necessary messages from part table''' 
    
    p_table = p_table.sort_values(by = "NORM_DSCORE", ascending = False)
    p_table["REP_SCORE"] = p_table["NORM_DSCORE"].values
    ptable_nrts = p_table["NORM_RT"].values
    ptable_fids = p_table["ID"].values
    ptable_rids = p_table["RUN_ID"].values
    ptable_rss = p_table["REP_SCORE"].values
    ptable_fid2repss = p_table.set_index('ID')['REP_SCORE'].to_dict()
    ptable_fid2nrts = p_table.set_index('ID')['NORM_RT'].to_dict()

    pt_bm = PT_BM(ptable_fid2repss, ptable_fid2nrts)
    
    p_nrt_list = np.arange(ptable_nrts.min() - 0.1, ptable_nrts.max() + 0.1, nrt_intervel)

    # 根据 NRT 距离以及REP_SCORE 计算mr score
    for nrt in p_nrt_list:
        nrt_label = np.abs(ptable_nrts - nrt) < nrt_width
        if np.sum(nrt_label) <= 1:
            pt_bm.mrss_list.append(0)
            pt_bm.mrs_list.append([])
            pt_bm.rid_list.append([])
            pt_bm.fid_list.append([])
            pt_bm.rid_u_list.append([])
            pt_bm.fid_u_list.append([])
            continue

        # ptp_mrs = rbf_kernel([[nrt]], ptable_nrts[nrt_label][:, None].tolist(), 1 / (2 * nrt_width))[0] * ptable_rss[nrt_label]
        ptp_mrs = rbf_kernel_numba_matrix(np.array([nrt]), ptable_nrts[nrt_label], nrt_width / 1.6)[0] * ptable_rss[nrt_label]
        ptp_rids = ptable_rids[nrt_label]
        ptp_fids = ptable_fids[nrt_label]

        # numba 加速
        ptp_u_inds = []
        for rid in np.unique(ptp_rids):
            ptp_rid_index = np.where(ptp_rids == rid)[0]
            if ptp_rid_index.shape[0] == 1:
                ptp_u_inds.append(ptp_rid_index[0])
            else:
                ptp_max_ind = np.argmax(ptp_mrs[ptp_rid_index])
                ptp_u_inds.append(ptp_rid_index[ptp_max_ind])
        pt_bm.fid_list.append(ptp_fids)
        pt_bm.rid_list.append(ptp_rids)           
        pt_bm.fid_u_list.append(ptp_fids[ptp_u_inds])
        pt_bm.rid_u_list.append(ptp_rids[ptp_u_inds])

        ptp_u_mrs = ptp_mrs[ptp_u_inds]

        pt_bm.mrs_list.append(ptp_u_mrs)
        pt_bm.mrss_list.append(ptp_u_mrs.sum())
    return pt_bm

def topN_mrgoup(mrss_list: list, n_mrg: int):

    # 确定前 n_mrg个候选mr Group 范围
    ptp_peaks, ptp_intens = find_peaks(mrss_list, distance = 8, height = 0.1)
    ptp_intens = ptp_intens["peak_heights"]
    if ptp_peaks.shape[0] == 0:
        return ptp_peaks, ptp_intens

    ptp_peak_rank = ptp_intens.argsort()[::-1][:n_mrg]
    ptp_peaks = ptp_peaks[ptp_peak_rank]
    ptp_intens = ptp_intens[ptp_peak_rank]
    return ptp_peaks, ptp_intens

class CMRG():
    def __init__(self, precur_id, rank, fids, rids, u_fids, u_rids):
        self.rank = rank
        self.fids = fids
        self.rids = rids
        self.init_u_fids = u_fids
        self.u_rids = u_rids
        self.nuf = len(u_fids)
        self.precur_id = precur_id
        
#     def set_init_messages(self, init_sub_js, sm_matrix, nrtdiff_matrix, reps_matrix, init_repss):
#         self.init_sub_js = init_sub_js
#         self.init_js = init_sub_js.sum()
        
#         self.init_sm_matrix = sm_matrix
#         self.init_nrtdiff_matrix = nrtdiff_matrix
#         self.init_reps_matrix = reps_matrix
#         self.init_repss = init_repss

        
    def set_iter_messages(self, iter_u_fids, iter_sub_js, sm_matrix, nrtdiff_matrix, reps_matrix, iter_repss):
        self.iter_u_fids = iter_u_fids
        self.iter_sub_js = iter_sub_js
        self.iter_js = iter_sub_js.sum()
        
        self.iter_sm_matrix = sm_matrix
        self.iter_nrtdiff_matrix = nrtdiff_matrix
        self.iter_reps_matrix = reps_matrix
        self.iter_repss = iter_repss
        

def return_cmrg(pt_bm, p_fid2xics, precur_id, min_nuf, nrt_width, debug_mode):
    
    cmrg_list = []
    
    # 比较每个mr group中peak group的最佳选择, 初始值为unique中的组合
    for n_ind, ptp_ind in enumerate(pt_bm.ptp_peaks):

        cmrg = CMRG(precur_id,
                    n_ind + 1,
                    pt_bm.fid_list[ptp_ind], 
                    pt_bm.rid_list[ptp_ind], 
                    pt_bm.fid_u_list[ptp_ind], 
                    pt_bm.rid_u_list[ptp_ind])

        if cmrg.nuf < min_nuf:
            continue

        u_fids_set = set(cmrg.init_u_fids)
        r_fids_label = [fid not in u_fids_set for fid in cmrg.fids]

        r_fids = cmrg.fids[r_fids_label]
        r_rids = cmrg.rids[r_fids_label]

        # initial
        initial_u_fids = cmrg.init_u_fids

        initial_mr_xics = np.array(itemgetter(*list(cmrg.init_u_fids))(p_fid2xics))
        initial_repss = np.array(itemgetter(*list(cmrg.init_u_fids))(pt_bm.fid2repss))
        initial_nrts = np.array(itemgetter(*list(cmrg.init_u_fids))(pt_bm.fid2nrts))

        sm_matrix = return_sm_matrix(initial_mr_xics)
        nrtdiff_matrix = rbf_kernel_numba_matrix(initial_nrts, initial_nrts, nrt_width / 1.6)
        reps_matrix = (initial_repss[:, np.newaxis] * initial_repss[np.newaxis, :])

        if debug_mode:
            assert sm_matrix.min() >= 0, f'min(sm_matrix) cannot be less than zero.'
            assert sm_matrix.max() <= 1, f'max(sm_matrix) cannot be more than one.'
            assert sm_matrix.shape[0] == sm_matrix.shape[1], f'sm matrix must have the same shape.'
            assert nrtdiff_matrix.min() >= 0, f'min(nrtdiff_matrix) cannot be less than zero.'
            assert nrtdiff_matrix.max() <= 1, f'max(nrtdiff_matrix) cannot be more than one.'
            assert nrtdiff_matrix.shape[0] == nrtdiff_matrix.shape[1], f'nrtdiff matrix must have the same shape.'
            assert reps_matrix.min() >= 0, f'min(reps_matrix) cannot be less than zero.'
            assert reps_matrix.shape[0] == reps_matrix.shape[1], f'reps matrix must have the same shape.'

        initial_sub_js = (sm_matrix * nrtdiff_matrix * reps_matrix).sum(axis = 0)
        # cmrg.set_init_messages(initial_sub_js, sm_matrix, nrtdiff_matrix, reps_matrix, initial_repss)

        # iteration
        iter_mr_xics = np.array(initial_mr_xics.tolist())
        iter_repss = np.array(initial_repss.tolist())
        iter_nrts = np.array(initial_nrts.tolist())
        iter_sub_js = np.array(initial_sub_js.tolist())
        iter_u_fids = np.array(cmrg.init_u_fids.tolist())

        for r_ind in range(len(r_fids)):
            rfid = r_fids[r_ind]
            rrid = r_rids[r_ind]
            rxic = p_fid2xics[rfid]
            rnrt = pt_bm.fid2nrts[rfid]
            rreps = pt_bm.fid2repss[rfid]

            r_location = np.where(cmrg.u_rids == rrid)[0][0]
            r_sm_vector = return_sm_vector(iter_mr_xics, rxic, r_location)
            r_nrtdiff_vector = rbf_kernel_numba_vector(iter_nrts, rnrt, r_location, nrt_width / 1.6)
            r_reps_vector = return_reps_vector(iter_repss, rreps, r_location)
            r_js = (r_sm_vector * r_nrtdiff_vector * r_reps_vector).sum()


            if r_js > iter_sub_js[r_location]:
                iter_mr_xics[r_location] = rxic
                iter_repss[r_location] = rreps
                iter_nrts[r_location] = rnrt

                iter_sm_matrix = return_sm_matrix(iter_mr_xics)
                iter_nrtdiff_matrix = rbf_kernel_numba_matrix(iter_nrts, iter_nrts, nrt_width / 1.6)
                iter_reps_matrix = (iter_repss[:, np.newaxis] * iter_repss[np.newaxis, :])

                iter_sub_js = (iter_sm_matrix * iter_nrtdiff_matrix * iter_reps_matrix).sum(axis = 0)

        cmrg.set_iter_messages(iter_u_fids, iter_sub_js, sm_matrix, nrtdiff_matrix, reps_matrix, iter_repss)
        cmrg_list.append(cmrg)
    return cmrg_list


def get_ndscore(f2nds_cur, feature_id) -> float:
    
    # with f2nds_conn.begin(write=False) as txn:
    #     value = txn.get(str(feature_id).encode('utf-8'))
    #     return(float(value.decode('utf-8')))

    value = f2nds_cur.get(str(feature_id).encode('utf-8'))
    return(float(value.decode('utf-8')))

def get_rid_native2chromid_db(rid_native2chromid_fpath: dict) -> dict:
    
    rid_native2chromid_db = {}
    rid_native2chromid_cur = {}   
    for k, v in rid_native2chromid_fpath.items():
        rid_native2chromid_db[k] = lmdb.open(v, create=False)
        rid_native2chromid_cur[k] = rid_native2chromid_db[k].begin(write=False)
    return rid_native2chromid_db, rid_native2chromid_cur


def get_cmrg_messages(precur_ids: list, db_fpath: str, feature2ndscore_fpath: str, rid_native2chromid_fpath: dict,
                      pr2tr_id_map: dict, rid2rn: dict, rn2chrom_fpath: dict, nrt_intervel: float, nrt_width: float,
                      n_mrg: int, min_nuf: int, logger_n, debug_mode, results_collector, counter, num_precursor, logger):
    # m_conn = sqlite3.connect(db_fpath)
    m_conn = sqlite3.connect(f'file:{db_fpath}?mode=ro', uri=True)
    # m_conn.close()
    results = []
    f2nds_conn = lmdb.open(feature2ndscore_fpath)
    f2nds_cur = f2nds_conn.begin(write=False)
    rid_native2chromid_db, rid_native2chromid_cur = get_rid_native2chromid_db(rid_native2chromid_fpath)
    rid2chrom_conn = get_rid2chrom_conn(rid2rn, rn2chrom_fpath)
    
    for prec_id in precur_ids:
        counter.value += 1
        if counter.value % logger_n == 0:
            logger.info(f"( {counter.value} / {num_precursor}) precursor has Calculated...")
        
        p_table = pd.read_sql_query(f'SELECT * FROM FEATURE WHERE PRECURSOR_ID = {prec_id}', m_conn)
        if p_table.shape[0] == 0:
            continue
        ndscores = []
        for f_id in p_table["ID"]:
            ndscores.append(get_ndscore(f2nds_cur, f_id))
        p_table["NORM_DSCORE"] = ndscores

        p_u_rids = p_table["RUN_ID"].unique()
        p_prec_ids = p_table["PRECURSOR_ID"].unique()
        if len(p_prec_ids) != 1:
            logger.error(f'Error.')
            raise
        p_prec_id = p_prec_ids[0]
        p_native_ids = ["%s_Precursor_i0"%p_prec_id] + pr2tr_id_map[prec_id]
        p_rid2xics, p_rid2rts = organize_chroms(p_u_rids, rid2chrom_conn, p_native_ids, rid_native2chromid_cur, logger)
        p_fid2xics, p_xic_unified_length = return_xics_by_fid(p_table, p_rid2xics, p_rid2rts, logger)

        # p_table, nrt_intervel, nrt_width
        pt_bm = return_bm_from_pt(p_table, nrt_intervel, nrt_width)
        # n_mrg
        pt_bm.ptp_peaks, pt_bm.ptp_intens = topN_mrgoup(pt_bm.mrss_list, n_mrg)
        if pt_bm.ptp_peaks.shape[0] == 0:
            continue
        cmrg_list = return_cmrg(pt_bm, p_fid2xics, prec_id, min_nuf, nrt_width, debug_mode)
        results.extend(cmrg_list)
    results_collector.extend(results)
    f2nds_cur.abort()
    f2nds_conn.close()
    m_conn.close()
    for k, v in rid_native2chromid_cur.items():
        v.abort()
    for k, v in rid_native2chromid_db.items():
        v.close()