import sqlite3
import pandas as pd

def get_target_os_features(table_name, suffix, target_fids, db_fpath):
    conn = sqlite3.connect(db_fpath)

    cur = conn.cursor()   
    # cur.execute(f'CREATE INDEX IF NOT EXISTS idx_feature_id ON {table_name} (FEATURE_ID)')
    cur.execute(f"PRAGMA table_info({table_name})")
    columns_info = cur.fetchall()
    cur.close()


    # 提取以'VAR'开头的列名
    var_columns = [info[1] for info in columns_info if info[1].startswith('VAR')]

    ex_columns = ["FEATURE_ID"] + var_columns
    ex_columns_str = ', '.join(ex_columns)
    tb_columns = ["FEATURE_ID"] + [col + suffix for col in var_columns]

    feature_pd = pd.read_sql(f'SELECT {ex_columns_str} FROM {table_name}', conn)
    conn.close()
    feature_pd.columns = tb_columns
    feature_pd = feature_pd.set_index("FEATURE_ID", drop = True)
    feature_pd = feature_pd.loc[target_fids, :]
    return feature_pd

def get_os_features(target_fids, db_fpath):
    
    ms1_feature_pd = get_target_os_features("FEATURE_MS1", "_MS1", target_fids, db_fpath)
    ms2_feature_pd = get_target_os_features("FEATURE_MS2", "_MS2", target_fids, db_fpath)
    feature_pd = pd.concat([ms1_feature_pd, ms2_feature_pd], axis=1)
    feature_pd = feature_pd.dropna(axis='columns', how='all')
    return feature_pd