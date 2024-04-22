import sqlite3
import numpy as np
import pandas as pd


def return_mr_features(results):
    mr_iter_messages = []
    for res in results:
        mr_iter_messages.append(res)
      
    columns = ["PRECURSOR_ID", "FEATURE_ID"] + ["MR_SCORE_%s"%(i+1) for i in range(18)]

    mr_iter_features = pd.DataFrame(mr_iter_messages, columns = columns)

    return mr_iter_features

def return_precid2decoy(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    precid2decoy = pd.read_sql(f'SELECT ID, DECOY FROM PRECURSOR', conn)
    conn.close()
    precid2decoy = precid2decoy.set_index("ID")["DECOY"].to_dict()
    data["decoy"] = data["PRECURSOR_ID"].apply(lambda x : precid2decoy[x])
    return data

def return_fid2runid(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    fid2runid = pd.read_sql(f'SELECT ID, RUN_ID FROM FEATURE', conn)
    conn.close()
    fid2runid = fid2runid.set_index("ID")["RUN_ID"].to_dict()
    data["run_id"] = data["FEATURE_ID"].apply(lambda x : fid2runid[x])
    return data

def return_fid2rt(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    fid2rt = pd.read_sql(f'SELECT ID, EXP_RT FROM FEATURE', conn)
    conn.close()
    fid2rt = fid2rt.set_index("ID")["EXP_RT"].to_dict()
    data["RT"] = data["FEATURE_ID"].apply(lambda x : fid2rt[x])
    return data

def return_fid2drt(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    fid2drt = pd.read_sql(f'SELECT ID, DELTA_RT FROM FEATURE', conn)
    conn.close()
    fid2drt = fid2drt.set_index("ID")["DELTA_RT"].to_dict()
    data["delta_rt"] = data["FEATURE_ID"].apply(lambda x : fid2drt[x])
    return data

def return_fid2leftwidth(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    fid2leftwidth = pd.read_sql(f'SELECT ID, LEFT_WIDTH FROM FEATURE', conn)
    conn.close()
    fid2leftwidth = fid2leftwidth.set_index("ID")["LEFT_WIDTH"].to_dict()
    data["leftWidth"] = data["FEATURE_ID"].apply(lambda x : fid2leftwidth[x])
    return data

def return_fid2rightwidth(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    fid2rightwidth = pd.read_sql(f'SELECT ID, RIGHT_WIDTH FROM FEATURE', conn)
    conn.close()
    fid2rightwidth = fid2rightwidth.set_index("ID")["RIGHT_WIDTH"].to_dict()
    data["rightWidth"] = data["FEATURE_ID"].apply(lambda x : fid2rightwidth[x])
    return data

def return_fid2intensity(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    fid2intensity = pd.read_sql(f'SELECT FEATURE_ID, AREA_INTENSITY FROM FEATURE_MS2', conn)
    conn.close()
    fid2intensity = fid2intensity.set_index("FEATURE_ID")["AREA_INTENSITY"].to_dict()
    data["Intensity"] = data["FEATURE_ID"].apply(lambda x : fid2intensity[x])
    return data

def return_fid2appar(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    fid2appar = pd.read_sql(f'SELECT FEATURE_ID, AREA_INTENSITY FROM FEATURE_MS1', conn)
    conn.close()
    fid2appar = fid2appar.set_index("FEATURE_ID")["AREA_INTENSITY"].to_dict()
    data["aggr_prec_Peak_Area"] = data["FEATURE_ID"].apply(lambda x : fid2appar[x])
    return data

def return_fid2appap(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    fid2appap = pd.read_sql(f'SELECT FEATURE_ID, APEX_INTENSITY FROM FEATURE_MS1', conn)
    conn.close()
    fid2appap = fid2appap.set_index("FEATURE_ID")["APEX_INTENSITY"].to_dict()
    data["aggr_prec_Peak_Apex"] = data["FEATURE_ID"].apply(lambda x : fid2appap[x])
    return data

def initial_format(db_fpath, mr_iter_features, os_feature_pd):

    mr_iter_features = return_precid2decoy(db_fpath, mr_iter_features)
    mr_iter_features = return_fid2runid(db_fpath, mr_iter_features)
    mr_iter_features = return_fid2rt(db_fpath, mr_iter_features)
    mr_iter_features = return_fid2drt(db_fpath, mr_iter_features)
    mr_iter_features = return_fid2leftwidth(db_fpath, mr_iter_features)
    mr_iter_features = return_fid2rightwidth(db_fpath, mr_iter_features)
    mr_iter_features = return_fid2intensity(db_fpath, mr_iter_features)
    mr_iter_features = return_fid2appar(db_fpath, mr_iter_features)
    mr_iter_features = return_fid2appap(db_fpath, mr_iter_features)

    mr_iter_features = mr_iter_features.set_index('FEATURE_ID')
    mr_iter_features = pd.concat([mr_iter_features, os_feature_pd], axis=1)
    return mr_iter_features


def return_precid2pepid(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    precid2pepid = pd.read_sql(f'SELECT PRECURSOR_ID, PEPTIDE_ID FROM PRECURSOR_PEPTIDE_MAPPING', conn)
    conn.close()
    precid2pepid = precid2pepid.set_index("PRECURSOR_ID")["PEPTIDE_ID"].to_dict()
    data["PEPTIDE_ID"] = data["PRECURSOR_ID"].apply(lambda x : precid2pepid[x])
    return data

def return_pepid2protid(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    pepid2protid = pd.read_sql(f'SELECT PEPTIDE_ID, PROTEIN_ID FROM PEPTIDE_PROTEIN_MAPPING', conn)
    conn.close()
    pepid2protid = pepid2protid.set_index("PEPTIDE_ID")["PROTEIN_ID"].to_dict()
    data["PROTEIN_ID"] = data["PEPTIDE_ID"].apply(lambda x : pepid2protid[x])
    return data

def return_protid2prot(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    protid2prot = pd.read_sql(f'SELECT ID, PROTEIN_ACCESSION FROM PROTEIN', conn)
    conn.close()
    protid2prot = protid2prot.set_index("ID")["PROTEIN_ACCESSION"].to_dict()
    data["ProteinName"]  = data["PROTEIN_ID"].apply(lambda x : protid2prot[x])
    return data

def return_precid2pmz(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    precid2pmz = pd.read_sql(f'SELECT ID, PRECURSOR_MZ FROM PRECURSOR', conn)
    conn.close()
    precid2pmz = precid2pmz.set_index("ID")["PRECURSOR_MZ"].to_dict()
    data["mz"] = data["PRECURSOR_ID"].apply(lambda x : precid2pmz[x])
    return data

def return_precid2pcharge(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    precid2pcharge = pd.read_sql(f'SELECT ID, CHARGE FROM PRECURSOR', conn)
    conn.close()
    precid2pcharge = precid2pcharge.set_index("ID")["CHARGE"].to_dict()
    data["Charge"] = data["PRECURSOR_ID"].apply(lambda x : precid2pcharge[x])
    return data

def return_runid2filename(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    runid2filename = pd.read_sql(f'SELECT ID, FILENAME FROM RUN', conn)
    conn.close()
    runid2filename = runid2filename.set_index("ID")["FILENAME"].to_dict()
    data["filename"] = data["run_id"].apply(lambda x : runid2filename[x])
    return data

def return_featid2irt(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    featid2irt = pd.read_sql(f'SELECT ID, NORM_RT FROM FEATURE', conn)
    conn.close()
    featid2irt = featid2irt.set_index("ID")["NORM_RT"].to_dict()
    data["iRT"] = data["FEATURE_ID"].apply(lambda x : featid2irt[x])
    return data

def return_pepid2seq(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    pepid2seq = pd.read_sql(f'SELECT ID, UNMODIFIED_SEQUENCE FROM PEPTIDE', conn)
    conn.close()
    pepid2seq = pepid2seq.set_index("ID")["UNMODIFIED_SEQUENCE"].to_dict()
    data["Sequence"] = data["PEPTIDE_ID"].apply(lambda x : pepid2seq[x])
    return data

def return_pepid2mseq(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    pepid2mseq = pd.read_sql(f'SELECT ID, MODIFIED_SEQUENCE FROM PEPTIDE', conn)
    conn.close()
    pepid2mseq = pepid2mseq.set_index("ID")["MODIFIED_SEQUENCE"].to_dict()
    data["FullPeptideName"] = data["PEPTIDE_ID"].apply(lambda x : pepid2mseq[x])
    return data

def return_precid2assayirt(db_fpath, data):
    conn = sqlite3.connect(db_fpath)
    precid2assayirt = pd.read_sql(f'SELECT ID, LIBRARY_RT FROM PRECURSOR', conn)
    conn.close()
    precid2assayirt = precid2assayirt.set_index("ID")["LIBRARY_RT"].to_dict()
    data["assay_iRT"] = data["PRECURSOR_ID"].apply(lambda x : precid2assayirt[x])
    return data


def output_format(db_fpath, mr_iter_res):
    mr_iter_res = mr_iter_res.reset_index()
    mr_iter_res = return_precid2pepid(db_fpath, mr_iter_res)
    mr_iter_res = return_pepid2protid(db_fpath, mr_iter_res)
    mr_iter_res = return_protid2prot(db_fpath, mr_iter_res)
    # mr_iter_res = return_precid2grouplabel(db_fpath, mr_iter_res)
    mr_iter_res = return_precid2pmz(db_fpath, mr_iter_res)
    mr_iter_res = return_precid2pcharge(db_fpath, mr_iter_res)
    mr_iter_res = return_runid2filename(db_fpath, mr_iter_res)
    mr_iter_res = return_featid2irt(db_fpath, mr_iter_res)
    mr_iter_res = return_pepid2seq(db_fpath, mr_iter_res)
    mr_iter_res = return_pepid2mseq(db_fpath, mr_iter_res)
    mr_iter_res["assay_rt"] = mr_iter_res["RT"] - mr_iter_res["delta_rt"]

    mr_iter_res["id"] = mr_iter_res["FEATURE_ID"]
    mr_iter_res = return_precid2assayirt(db_fpath, mr_iter_res)
    mr_iter_res["delta_iRT"] = mr_iter_res["iRT"] - mr_iter_res["assay_iRT"]
    
    scored_columns = []
    bases_columns = []
    for col in mr_iter_res.columns:
        if col.startswith("VAR_") or col.startswith("MR_"):
            scored_columns.append(col)
        else:
            bases_columns.append(col)

    mr_iter_res = mr_iter_res.loc[:, bases_columns]
    return mr_iter_res