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

# @file  
# @brief 
# @author Ryan Feehan <RFeehan93@gmail.com>

import json
import pandas as pd
import MAHOMES_II as MHMII
import sys

job_dir = str(sys.argv[1])

features = pd.read_csv("%s/features.csv"%(job_dir))

## place bad sites to the side while other sites get predictions
no_preds = features.loc[features['bad_site']==True].copy()
if 0<len(no_preds): 
    no_preds = no_preds[['struc_id', 'note']].copy()
    no_preds.rename(columns={"struc_id":"jobid", "note":"prediction"}, inplace=True)
else:
    no_preds=pd.DataFrame(columns=["jobid", "prediction"])


## make predictions or empty df if no sites had successful feature calc.
features = features.loc[features['bad_site']==False].copy()
if 0<len(features):
    ## set all site info to the side so that it can be added to the returned predictions from MAHOMES II
    save_preds = features.copy()
    save_preds['metal1_seqID']=save_preds['metal1_seqID'].astype(int)
    save_preds.set_index(['struc_id', 'metal1_resName', 'metal1_seqID', 'metal1_serial'], inplace=True)
    save_preds[['junk', 'site#']]=save_preds['SITE_ID'].str.rsplit("_", n=1, expand=True)
    
    for col in ['metal2_resName', 'metal2_seqID', 'metal3_resName', 'metal3_seqID', 'metal4_resName', 'metal4_seqID']:
        if col not in save_preds.columns.tolist(): save_preds[col]=""
        elif col in ['metal2_seqID', 'metal3_seqID', 'metal4_seqID']: save_preds[col]=save_preds[col].astype(int)

    save_preds = save_preds[['site#', 'metal2_resName', 'metal2_seqID', 'metal3_resName', 'metal3_seqID', 'metal4_resName', 'metal4_seqID']].copy()
    # make predictions    
    preds = MHMII.make_predictions_with_saved_MAHOMES_II(features)
    save_preds = save_preds.merge(preds[['prediction', 'bool_pred']].copy(), left_index=True, right_index=True)

    ## rename as needed
    save_preds.reset_index(inplace=True)
    save_preds.rename(columns={"struc_id":"jobid", "prediction":"percent catalytic predictions"}, inplace=True)
    for i in range(1,5):
        save_preds.rename(columns={"metal%s_resName"%(str(i)):"Name%s"%(str(i)), "metal%s_seqID"%(str(i)):"Res#%s"%(str(i))}, inplace=True)
    ## turn prediction value into words    
    save_preds['prediction']="Not Catalytic"
    save_preds.loc[save_preds["bool_pred"]==True,'prediction']="Catalytic"
    save_preds["percent catalytic predictions"]=save_preds["percent catalytic predictions"]*100

else:
    save_preds = pd.DataFrame(columns=['jobid','input file', 'site#'
        ,'Name1', 'Res#1', 'Name2', 'Res#2','Name3', 'Res#3', 'Name4', 'Res#4' 
        ,"percent catalytic predictions",'prediction'])

## add back input structures that failed along the way
save_preds = pd.concat([save_preds, no_preds], ignore_index=True, join='outer')
save_preds.fillna("", inplace=True)

## add user input file name
with open("%s/input-legend.json"%(job_dir), "r") as f:
    name_dict = json.loads(f.read())
save_preds['jobid']=save_preds['jobid']+".pdb"
save_preds['input file']=save_preds['jobid'].map(name_dict)
save_preds['input file'] = save_preds['input file'].str[:-4]

## reorder, which will change when I change col names!
save_preds = save_preds[[
    'jobid','input file', 'site#'
    ,'Name1', 'Res#1', 'Name2', 'Res#2' #,'metal1_resName', 'metal1_seqID', 'metal2_resName', 'metal2_seqID'
    ,'Name3', 'Res#3', 'Name4', 'Res#4' #,'metal3_resName', 'metal3_seqID', 'metal4_resName', 'metal4_seqID'
    ,"percent catalytic predictions",'prediction']].copy()
for col in ["Res#1", "Res#2", "Res#3", "Res#4"]:
    save_preds[col]=pd.to_numeric(save_preds[col],  errors='coerce', downcast="integer")

save_preds.to_csv("%s/predictions.csv"%(job_dir), index=False)

