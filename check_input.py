import sqlite3

def check_table_exists(db_conn, table_name: str) -> bool:
    """
    Check if a table exists in the SQLite database.
    :param db_conn: a sqlite3.Connection object
    :param table_name: name of the table to check
    :return: True if the table exists, False otherwise
    """
    # 创建一个新的cursor对象
    cur = db_conn.cursor()

    # 查询sqlite_master表检查特定表名是否存在
    cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?", (table_name,))

    # 如果返回结果不为空，说明表存在
    if cur.fetchone():
        cur.close()
        return True

    cur.close()
    return False

def check_column_exists(db_conn, table_name: str, column_name: str):
    """
    Check if a column exists in a specific table of the SQLite database.
    :param db_conn: a sqlite3.Connection object
    :param table_name: name of the table to check
    :param column_name: name of the column to check for
    :return: True if the column exists, False otherwise
    """
    cur = db_conn.cursor()

    # 查询特定表的列信息
    cur.execute(f"PRAGMA table_info({table_name})")

    # 获取查询结果中的所有列名
    columns = [info[1] for info in cur.fetchall()]
    cur.close()

    return column_name in columns

def check_db(db_file: str, logger) -> None:
    
    checked_m = {}
    checked_m["SCORE_MS2"] = ["FEATURE_ID", "SCORE"]
    checked_m["FEATURE"] = ["ID", "RUN_ID", "PRECURSOR_ID", "EXP_RT", "NORM_RT", "LEFT_WIDTH", 
                            "RIGHT_WIDTH"]
    checked_m["PRECURSOR"] = ["ID", "GROUP_LABEL", "DECOY"]
    checked_m["TRANSITION_PRECURSOR_MAPPING"] = ["TRANSITION_ID", "PRECURSOR_ID"]

    conn = sqlite3.connect(db_file)

    for table, cols in checked_m.items():
        if not check_table_exists(conn, table):
            logger.error(f'An error occurred: {table}(table) does not exist in {db_file}')
            conn.close()
            raise
        for col in cols:
            if not check_column_exists(conn, table, col):
                logger.error(f'An error occurred: {col}(column) does not exist in {table}(table)')
                conn.close()
                raise
    conn.close()