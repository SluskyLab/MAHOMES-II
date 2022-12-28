#    MAHOMES II app
#    Copyright (C) 2021 University of Kansas
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published
#    by the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# @file   FeatureCalculations/CalcFeatureDROPP.py
# @brief  Calculates the DROPP for input feature values of enzyme and non-enzyme sites


import pandas as pd
import numpy as np
from scipy import stats
from scipy import integrate
import math
from collections import Counter
import os
import sys 


ACPT_STDS = 10

not_features = ["Enzyme", "Set", "SITE_ID", "struc_id", "bad","bad_site", "note"]
for x in range(1,5):
    not_features.append("metal%s_seqID"%(str(x)))
    not_features.append("metal%s_serial"%(str(x)))
    not_features.append("metal%s_resName"%(str(x)))
    
    
def calculate_DROPP(features_df, acceptable_stds=ACPT_STDS, output_file=None):
    sim_values = []
    for term in features_df.columns.tolist():
        numeric_feature=True
        try:
            features_df[term] = features_df[term].astype(float)
        except:
            numeric_feature=False
        if numeric_feature:
            if term not in not_features:
                df = features_df[['Enzyme', term]].copy()
                ## no variation means all values are the same and we can just set sim.=1
                std = df[term].std()
                if 0==std:
                    jac_score=1; auc=1; num_unique_vals=1; removed_outliers=0
                else:
                    removed_outliers = 0
                    unique_vals = df[term].unique().tolist()
                    if len(unique_vals)>2:
                        removed_outliers = len(df[np.abs(df[term]-df[term].mean()) > (acceptable_stds*std)])
                        df = df[np.abs(df[term]-df[term].mean()) <= (acceptable_stds*std)]
                    unique_vals = df[term].unique().tolist()
                    num_unique_vals = len(unique_vals)

                    totalA = len(df[df["Enzyme"] == True])
                    totalB = len(df[df["Enzyme"] == False])
                    this_prob = Counter(df[term].loc[df["Enzyme"] == True])
                    bg_prob = Counter(df[term].loc[df["Enzyme"] == False])
                    shared_keys = sorted(set(this_prob.keys()).union(bg_prob.keys()))
                    counts = np.zeros([2,len(shared_keys)])
                    for key in range(0,len(shared_keys)):
                        counts[0,key] = this_prob[shared_keys[key]]
                        counts[1,key] = bg_prob[shared_keys[key]]
                    counts /= counts.sum(axis=1, keepdims=True)
                    jac_score = np.sum(counts.min(axis=0)) / np.sum(counts.max(axis=0))
                    try:
                        kde1 = stats.gaussian_kde(df[term].loc[df["Enzyme"] == True].dropna())
                        kde2 = stats.gaussian_kde(df[term].loc[df["Enzyme"] == False].dropna())
                        this_min = np.min(df[term])
                        this_max = np.max(df[term])
                        this_interval = np.linspace(this_min, this_max, 2**10+1)
                        cat_pts = kde1.evaluate(this_interval)
                        non_cat_pts = kde2.evaluate(this_interval)
                        int_pts = np.minimum(cat_pts, non_cat_pts)
                        step_size = this_interval[1] - this_interval[0]
                        if step_size<0:step_size*=-1
                        auc = integrate.romb(int_pts, dx=step_size)
                    except:
                        auc=1
                        print("%s fail auc calc."%term)
                sim_values.append([term,jac_score, auc, num_unique_vals, removed_outliers ])

    feat_signif = pd.DataFrame.from_records(sim_values, columns = ["Feature", 'Jaccard Similarity', "AUC Similarity", "Different Values", "Num. removed outliers"])

    feat_signif['DROPP']=feat_signif['Jaccard Similarity']
    feat_signif.loc[feat_signif['Different Values']>20, 'DROPP']=feat_signif['AUC Similarity']
    feat_signif = feat_signif.sort_values(['DROPP'],ascending=[True]).reset_index(drop=True)
    if None!=output_file:
        feat_signif.to_csv(output_file)
        
    return(feat_signif)



