import sqlite3
import lmdb
import pandas as pd
import numpy as np

from operator import itemgetter
# from scipy.signal import find_peaks
from scipy.spatial.distance import cosine

from database import get_rid2chrom_conn
# from utils import return_sm_matrix, return_sm_vector, rbf_kernel_numba_matrix, rbf_kernel_numba_vector, return_reps_vector
from utils import pearson_corrcoef_between_arrays
from chromatographic import organize_chroms, return_xics_by_fid


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
                      pr2tr_id_map: dict, rid2rn: dict, rn2chrom_fpath: dict, nrt_width: float,
                      logger_n, results_collector, counter, num_precursor, logger):

    m_conn = sqlite3.connect(f'file:{db_fpath}?mode=ro', uri=True)
    # m_conn.close()
    results = []
    f2nds_conn = lmdb.open(feature2ndscore_fpath)
    f2nds_cur = f2nds_conn.begin(write=False)
    rid_native2chromid_db, rid_native2chromid_cur = get_rid_native2chromid_db(rid_native2chromid_fpath)
    rid2chrom_conn = get_rid2chrom_conn(rid2rn, rn2chrom_fpath)

    results = []

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
            logger.error(f'Error')
            raise
        p_prec_id = p_prec_ids[0]
        p_native_ids = ["%s_Precursor_i0"%p_prec_id] + pr2tr_id_map[prec_id]
        p_rid2xics, p_rid2rts = organize_chroms(p_u_rids, rid2chrom_conn, p_native_ids, rid_native2chromid_cur, logger)
        p_fid2xics, p_xic_unified_length = return_xics_by_fid(p_table, p_rid2xics, p_rid2rts, logger)

        run_maxids = p_table.groupby('RUN_ID')['NORM_DSCORE'].idxmax()

        for p_best_peak_ind in run_maxids:
            # p_best_peak_ind = np.argmax(p_table["NORM_DSCORE"].values)
            m_nds = p_table.iloc[p_best_peak_ind]["NORM_DSCORE"]
            m_nrt = p_table.iloc[p_best_peak_ind]["NORM_RT"]
            m_id = p_table.iloc[p_best_peak_ind]["ID"]
            p_nrts_diff = np.abs(p_table["NORM_RT"].values - m_nrt)
            p_nrt_labels = p_nrts_diff < nrt_width / 2
            sp_table = p_table.loc[p_nrt_labels, :].copy()
            if sp_table.shape[0] == 0:
                continue

            m_xics = p_fid2xics[m_id][: 7]
            m_intens = m_xics.max(axis = 1)
            if (m_intens == 0).all():
                continue
            smss = []
            for fid in sp_table["ID"]:
                pear_matrix = pearson_corrcoef_between_arrays(m_xics[: 7], p_fid2xics[fid][: 7])
                pear_diag = np.diag(pear_matrix)
                smss.append(pear_diag.sum())
            sp_table["SMSS"] = smss
            sp_table = sp_table.sort_values(["RUN_ID", "SMSS"], ascending = [False, False])
            s_sp_table = sp_table.drop_duplicates("RUN_ID", keep = "first")
            if s_sp_table.shape[0] == 0:
                continue
            for ind in s_sp_table.index:
                fid, run_id, smss, nrt = s_sp_table.loc[ind, ["ID", "RUN_ID", "SMSS", "NORM_RT"]]
                if smss < 3.:
                    continue
                else:
                    c_xics = p_fid2xics[fid][: 7]
                    c_intens = c_xics.max(axis = 1)
                    if (c_intens == 0).all():
                        continue
                    pear_matrix = pearson_corrcoef_between_arrays(m_xics, c_xics)
                    tmp = np.zeros(17)
                    diag_vector = np.diag(pear_matrix)
                    for ii, pears in enumerate(pear_matrix):
                        tmp[ii] = pears[np.argsort(pears)[::-1][1]]
                    tmp[7: 7 + len(pear_matrix.mean(axis = 0))] = pear_matrix.mean(axis = 0)
                    tmp[14] = 1 - cosine(m_intens, c_intens)
                    tmp[15] = np.abs(m_nrt - nrt)
                    tmp[16] = m_nds
                    results.append([prec_id, fid, smss] + tmp.tolist())

    results_collector.extend(results)
    f2nds_cur.abort()
    f2nds_conn.close()
    m_conn.close()
    for k, v in rid_native2chromid_cur.items():
        v.abort()
    for k, v in rid_native2chromid_db.items():
        v.close()