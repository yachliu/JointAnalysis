import os
import sqlite3

import numpy as np
import pandas as pd

from tqdm import tqdm
from bisect import bisect
from loguru  import logger

from multiprocessing import Pool, cpu_count, Manager, Process, Value

# pd.set_option('display.max_columns', None)
from check_input import check_db
from preprocessing import build_feature2ndscore, get_db_rid2rn, get_db_rn2fpath, return_pr2tr_id_map, return_nrt_width
from database import get_rid2chrom_conn, close_rid2chrom_conn, get_run_native2chrom_fpath
from mrgroup import get_cmrg_messages
from format_data import return_mr_features, initial_format, output_format
from openswath_feature import get_os_features
from discriminate import calc_score_cut, calc_results
from reports import stats

def mrgd(db_fpath, chrom_dpath, work_dpath, n_threads, map_size, tfdr, nrt_width_percent, seed):
    map_size = 2 ** map_size
    if not os.path.exists(work_dpath):
        os.makedirs(work_dpath)
    log_fpath = os.path.join(work_dpath, "JointAnalysis.log")
    logger.add(log_fpath, format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}", mode="w")
    logger.info(f'JointAnalysis Workflow')


    logger.info(f'Check db_infile: {db_fpath}')
    check_db(db_fpath, logger)

    feature2ndscore_fpath = os.path.join(work_dpath, "feature2ndscore.db")
    logger.info(f'Save ndscores to db: {feature2ndscore_fpath}')
    build_feature2ndscore(db_fpath, feature2ndscore_fpath, map_size)

    logger.info(f'Organize the necessary inputs')
    rid2rn = get_db_rid2rn(db_fpath)
    rid_list = [k for k in rid2rn.keys()]
    rn2chrom_fpath = get_db_rn2fpath(chrom_dpath, "sqMass")
    pr2tr_id_map = return_pr2tr_id_map(db_fpath)

    nrt_width = return_nrt_width(db_fpath, nrt_width_percent)

    logger.info(f'Save nativeID2chromID')
    rid2chrom_conn = get_rid2chrom_conn(rid2rn, rn2chrom_fpath)
    rid_native2chromid_fpath = get_run_native2chrom_fpath(rid2chrom_conn, work_dpath, map_size / 8)
    close_rid2chrom_conn(rid2chrom_conn)

    logger.info(f'Get MRGroup')
    m_conn = sqlite3.connect(db_fpath)
    m_cur = m_conn.cursor()
    m_cur.execute(f'SELECT ID FROM PRECURSOR')

    precursor_ids = np.array(m_cur.fetchall()).squeeze()
    m_cur.close()
    m_conn.close()

    num_precursor = precursor_ids.shape[0]
    logger_n = 10 ** (len(str(num_precursor)) - 2)
    n_precur = num_precursor // n_threads
    precurs_list = [precursor_ids[i * n_precur : (i + 1) * n_precur].tolist() for i in range(n_threads)]
    _ = [precurs_list[i].append(precursor_ids[i + n_precur * n_threads]) for i in range(len(precursor_ids) - n_precur * n_threads)]

    results_collector = Manager().list()
    counter = Manager().Value('d',0)
    logger.info(f"( {counter.value} / {num_precursor}) precursor has Calculated...")
    extractors = []
    for precur_ids in precurs_list:
        p = Process(target = get_cmrg_messages, 
                    args =  (precur_ids, db_fpath, feature2ndscore_fpath, rid_native2chromid_fpath,
                             pr2tr_id_map, rid2rn, rn2chrom_fpath, nrt_width,
                             logger_n, results_collector, counter, num_precursor, logger, ))
        p.daemon = True
        extractors.append(p)
        p.start()
    for p in extractors:
        p.join()

    logger.info(f'Get MR features')
    mr_iter_features = return_mr_features(results_collector)
    del results_collector, counter

    logger.info(f'Get OS features')
    target_fids = mr_iter_features["FEATURE_ID"].values
    os_feature_pd = get_os_features(target_fids, db_fpath)

    logger.info(f'Initial format')
    mr_iter_features = initial_format(db_fpath, mr_iter_features, os_feature_pd)

    logger.info(f'Discriminate')
    ignored_columns = ["PRECURSOR_ID", "decoy", "run_id", "RT", "delta_rt", "rightWidth", "leftWidth", "Intensity", "aggr_prec_Peak_Area", "aggr_prec_Peak_Apex"]
    iter_mr_columns = [col for col in mr_iter_features.columns if col not in ignored_columns]

    mr_iter_res = calc_results(scored_columns = iter_mr_columns,
                               data_pd = mr_iter_features,
                               n_threads = n_threads,
                               seed = seed,
                               n_rawdata = len(rid2rn))

    logger.info(f'Output results')
    mr_iter_res = output_format(db_fpath, mr_iter_res)

    mr_iter_res = stats(mr_iter_res, "jd_score", logger)
    mr_iter_res = mr_iter_res.sort_values(["PRECURSOR_ID", "jd_score"], ascending=[True, False])
    mr_iter_res["peak_group_rank"] = mr_iter_res.groupby('PRECURSOR_ID').cumcount() + 1
    mr_iter_res["d_score"] = mr_iter_res["jd_score"]
    mr_iter_res["m_score"] = mr_iter_res["qvalue"]
    
    trans = []
    for _, (mseq, charge, decoy) in enumerate(zip(mr_iter_res["FullPeptideName"].values, mr_iter_res["Charge"].values, mr_iter_res["decoy"].values)):
        if decoy == 0:
            trans.append(mseq + "_" + str(charge))
        else:
            trans.append("DECOY_" + mseq + "_" + str(charge))
    mr_iter_res["transition_group_id"] = trans
    results_format = mr_iter_res.loc[:, ["transition_group_id",
                                        "decoy",
                                        "run_id",
                                        "filename", 
                                        "RT",
                                        "assay_rt",
                                        "delta_rt",
                                        "iRT",
                                        "assay_iRT",
                                        "delta_iRT",
                                        "id",
                                        "Sequence",
                                        "FullPeptideName",
                                        "Charge",
                                        "mz",
                                        "Intensity",
                                        "aggr_prec_Peak_Area",
                                        "aggr_prec_Peak_Apex",
                                        "leftWidth",
                                        "rightWidth",
                                        "peak_group_rank",
                                        "d_score",
                                        "m_score",
                                        "ProteinName"]]


    result_saved_fpath = os.path.join(work_dpath, "jointAnalysis_results.tsv")
    results_format.to_csv(result_saved_fpath, sep = "\t", index = False)
    logger.info(f'jointAnalysis_results: {result_saved_fpath}')

    logger.info(f'Done')





