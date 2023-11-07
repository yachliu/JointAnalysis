import os
import lmdb
import sqlite3
import pandas as pd


def get_rid2chrom_conn(rid2rn: dict, rn2fpath: dict) -> dict:
    rid2conn = {}
    for k,v in rid2rn.items():
        fpath = rn2fpath[v]
        # rid2conn[k] = sqlite3.connect(fpath)
        rid2conn[k] = sqlite3.connect(f'file:{fpath}?mode=ro', uri=True)
    return rid2conn

def close_rid2chrom_conn(rid2conn: dict) -> None:
    for k, v in rid2conn.items():
        v.close()

def save_run_native2chrom(saved_fpath, data_dict, map_size):
    
    n2c_db = lmdb.open(saved_fpath, create = True, map_size = map_size)
    with n2c_db.begin(write=True) as txn:
        for key, value in data_dict.items():
            txn.put(key.encode('utf-8'), value.encode('utf-8'))
    n2c_db.close()
    
def get_run_native2chrom_fpath(rid2chrom_conn, work_dpath: str, map_size) -> dict:
    
    rid_native2chromid_fpath = {}
    for rid, r_conn in rid2chrom_conn.items():
        c_table = pd.read_sql(f'SELECT NATIVE_ID, ID FROM CHROMATOGRAM', r_conn)
        c_table["ID"] = c_table["ID"].astype(str)
        c_table["NATIVE_ID"] = c_table["NATIVE_ID"].astype(str)
        run_native2chrom_id = dict(zip(c_table["NATIVE_ID"].values, c_table["ID"]))
        run_n2c_fpath = os.path.join(work_dpath, 'native2chromID_%s.db'%rid)
        save_run_native2chrom(run_n2c_fpath, run_native2chrom_id, map_size)
        rid_native2chromid_fpath[rid] = run_n2c_fpath
    return rid_native2chromid_fpath