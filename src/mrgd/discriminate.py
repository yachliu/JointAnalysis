import xgboost
import collections
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pulearn import BaggingPuClassifier
from xgboost.sklearn import XGBClassifier

def calc_score_cut(altsv, score_column, label_column, cut_off, smooth_factor = 0.01, plot = True, plot_name = None):  
    target = altsv[altsv[label_column] == 0]
    target = target.sort_values(by = score_column, ascending = False)
    decoy = altsv[altsv[label_column] != 0]
    decoy = decoy.sort_values(by = score_column, ascending = False)    
    def larger_count(alist, value):
        return list(map(lambda x : 1 if x >= value else 0, alist)).count(1)    
    target_scores = list(target[score_column])
    decoy_scores = list(decoy[score_column])
    fdrs = []
    n_fp = 0
    for i, cut in enumerate(target_scores):
        n_tp = i + 1
        while n_fp < len(decoy_scores) and decoy_scores[n_fp] >= cut:
            n_fp += 1
        fdr = n_fp / (n_fp + n_tp)
        fdrs.append(fdr)    
        
    cut_modified, times = 0, 0 
    final_cut = 0.7
    while cut_modified == 0 and times < 20 and int(smooth_factor * len(fdrs)) > 0:
        slider = collections.deque()
        for i in range(int(smooth_factor * len(fdrs))):
            slider.append(fdrs[i])            
        for i in range(int(smooth_factor * len(fdrs)), len(fdrs)):
            if slider[0] >= cut_off:
                if larger_count(slider, cut_off) >= len(slider) * 0.8:
                    final_cut = target_scores[i - 1]
                    cut_modified = 1
                    break
            slider.popleft()
            slider.append(fdrs[i])
        smooth_factor /= 1.2
        times += 1    
    if cut_modified == 0:
        print("   Warning: failed to calculate FDR. Used a heuristic cut-off value: %s" % final_cut)   
    if plot:
        plt.figure(figsize = (9.7, 7.4))
        sns.distplot(target[score_column], label = "target")
        sns.distplot(decoy[score_column], label = "decoy")
        plt.xticks(fontsize = 18)
        plt.yticks(fontsize = 18)
        plt.xlabel(score_column, fontsize = 22)
        plt.ylabel("Frequency", fontsize = 22)
        plt.legend(fontsize = 22)
        plt.axvline(x = final_cut, color = "red")
        plt.title("%.6f" % final_cut)
        plt.show()
    return fdrs, final_cut


def calc_results(scored_columns, data_pd, n_threads, seed, n_rawdata):

    max_depth = 5

    trains = data_pd.loc[:, scored_columns].values
    labels = data_pd["decoy"].values

    cl = XGBClassifier(n_jobs = n_threads, max_depth = max_depth, random_state = seed, gamma = 0.05, reg_alpha = 0.01,
                       reg_lambda = 0.1, subsample = 0.8 / n_rawdata, use_label_encoder=False, eval_metric="logloss")

    cl_pu = BaggingPuClassifier(cl)
    cl.fit(trains, labels)
    data_pd["jd_score"] = cl.predict_proba(data_pd.loc[:, scored_columns].values)[:, 0]
    
    scored_pd_iteration = data_pd.sort_values(by = ["run_id", "PRECURSOR_ID", "jd_score"], ascending = [True, True, False])
    scored_pd_iteration = scored_pd_iteration.drop_duplicates(["run_id", "PRECURSOR_ID"], keep = "first")

    # fdrs, final_cut = calc_score_cut(scored_pd_iteration, "JOINT_DS", "decoy", 0.005, smooth_factor = 0.01, plot = True)
    # results = scored_pd_iteration[scored_pd_iteration["JOINT_DS"] >= final_cut]
    # print(results['decoy'].value_counts())
    # results = results[results['decoy'] == 0]
    # print(results['species'].value_counts())
    return scored_pd_iteration