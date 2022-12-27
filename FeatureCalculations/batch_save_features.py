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
# @file   FeatureCalculations/batch_save_features.py 
# @brief  Runs and saves features for structure files in input directory.
# @author Ryan Feehan <RFeehan93@gmail.com>

import os
import glob
import pandas as pd
import get_pdb_features as features
import sys

## third party software has previously been ran for this job
job_dir = str(sys.argv[1]) 

def run_save_features(data_dir):
    with open("%s/batch_input.txt"%data_dir, 'r') as f:
        struc_list = f.read()
    print(struc_list)
    struc_list = struc_list.split("\n")
    print(struc_list)
    all_job_features = []
    for struc in struc_list:
        if (1<=len(struc)) and (os.path.isdir("%s/%s"%(data_dir, struc))): 
            struc_features = features.get_features_for_pdb(data_dir, struc) 
            if len(struc_features)>0:
                struc_features.set_index("SITE_ID", drop=True, inplace=True)
                all_job_features.append(struc_features)
            else:
                print("OH NO: failed %s without knowning why"%struc)
                ## TODO: add check here to see why, like were no metals present?
                ##       and return df w/ struc name and note the issue

    ## save calculations for all sites on dataset PDBs
    all_job_features = pd.concat(all_job_features, join='outer', axis=0)
    all_job_features.to_csv("%s/features.csv"%(data_dir))
    print("Finished features!", all_job_features.shape)
run_save_features(job_dir)

