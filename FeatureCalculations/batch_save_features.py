import os
import glob
import pandas as pd
import get_pdb_features as features
import sys

## third party software has previously been ran for this job
job_dir = str(sys.argv[1]) 

def run_save_features(data_dir):
    input_dirs = glob.glob("%s/*"%(data_dir))
    print(len(input_dirs))
    data_pdb_features = []
    for this_dir in input_dirs:
        if os.path.isdir(this_dir) == True: 
            pdb = this_dir.split("/")[-1]
            cur_pdb_features = features.get_features_for_pdb(data_dir, pdb, old_features=False)
            if len(cur_pdb_features)>0:
                cur_pdb_features.set_index("SITE_ID", drop=True, inplace=True)
                data_pdb_features.append(cur_pdb_features)
            else:
                print("OH NO: failed %s without knowning why %s"%pdb)

    ## save calculations for all sites on dataset PDBs
    data_pdb_features_df = pd.concat(data_pdb_features, join='outer', axis=0)
    data_pdb_features_df.to_csv("%s/features.csv"%(data_dir))
    print("Finished features!", data_pdb_features_df.shape)
run_save_features(job_dir)

