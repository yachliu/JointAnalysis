import os
import lmdb
import sqlite3
import pandas as pd
import numpy as np

from sklearn.preprocessing import StandardScaler
from scipy import stats

def build_feature2ndscore(db_fpath: str, saved_fpath: str, map_size) -> None:
    '''Feature ID 2 Normalize dscore'''

    conn = sqlite3.connect(db_fpath)
    s_table = pd.read_sql_query(f'SELECT FEATURE_ID, SCORE, RANK FROM SCORE_MS2', conn)
    scaler = StandardScaler()
    n_scores = scaler.fit_transform(s_table["SCORE"].values[:, np.newaxis]).squeeze().astype("float32")
    nn_scores = stats.norm.cdf(n_scores)
    features = s_table["FEATURE_ID"].values
    feature2ndscore = dict(zip(features, nn_scores))
    
    s_db = lmdb.open(saved_fpath, create = True, map_size = map_size)
    with s_db.begin(write=True) as txn:
        for key, value in feature2ndscore.items():
            txn.put(str(key).encode('utf-8'), str(value).encode('utf-8'))
    s_db.close()
    return None

def get_db_rid2rn(db_fpath: str) -> dict:
    m_conn = sqlite3.connect(db_fpath)
    rid2rn_pd = pd.read_sql('''SELECT * FROM RUN''', m_conn)
    m_conn.close()
    rid2rn = {}
    for row_id in rid2rn_pd.index:
        runid, fpath = rid2rn_pd.loc[row_id, :]
        base_name = os.path.basename(fpath).split(".")[0]
        rid2rn[runid] = base_name
    return rid2rn

def get_db_rn2fpath(directory, extension):
    """查找给定目录及其子目录下所有具有指定扩展名的文件并返回它们的绝对路径"""
    rn2fpath = {}
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(extension):
                absolute_path = os.path.abspath(os.path.join(root, file))
                rn2fpath[file.split(".")[0]] = absolute_path
    return rn2fpath

    
# def get_chromID_from_db(db, search_key):
#     # db = plyvel.DB(db_path, create_if_missing=False)
#     value = db.get(search_key.encode('utf-8'))
    
#     return value.decode('utf-8') if value else None

def return_pr2tr_id_map(db_fpath) -> dict:
    db_conn = sqlite3.connect(db_fpath)
    trans_prec_pd = pd.read_sql(f'SELECT * FROM TRANSITION_PRECURSOR_MAPPING', db_conn)
    p2t_id_map = trans_prec_pd.groupby('PRECURSOR_ID')['TRANSITION_ID'].apply(list).to_dict()
    db_conn.close()
    return p2t_id_map

def return_nrt_width(db_fpath, nrt_width_percent):
    db_conn = sqlite3.connect(db_fpath)
    norm_rts = pd.read_sql(f'SELECT NORM_RT FROM FEATURE', db_conn)["NORM_RT"].values
    nrt_width = (norm_rts.max() - norm_rts.min()) * nrt_width_percent
    db_conn.close()
    return nrt_width