import os
import glob
import pandas as pd
import get_pdb_features as features
import sys

## third party software has previously been ran for this job
job_dir = str(sys.argv[1]) 

def run_save_features(data_dir):
    with open("%s/batch_input.txt"%data_dir, 'r') as f:
        struc_list = f.readlines()
    all_job_features = []
    for struc in struc_list:
        struc_features = features.get_features_for_pdb(data_dir, struc)
        if len(struc_features)>0:
            struc_features.set_index("SITE_ID", drop=True, inplace=True)
            all_job_features.append(struc_features)
        else:
            print("OH NO: failed %s without knowning why %s"%struc)
            ## TODO: add check here to see why, like were no metals present?
            ##       and return df w/ struc name and note the issue

    ## save calculations for all sites on dataset PDBs
    all_job_features = pd.concat(all_job_features, join='outer', axis=0)
    all_job_features.to_csv("%s/features.csv"%(data_dir))
    print("Finished features!", all_job_features.shape)
run_save_features(job_dir)

