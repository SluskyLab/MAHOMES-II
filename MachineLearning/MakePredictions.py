import json
import pandas as pd
import MAHOMES_II as MHMII
import sys
import pandas as pd

job_dir = str(sys.argv[1])
features = pd.read_csv("%s/features.csv"%(job_dir))
preds = MHMII.make_predictions_with_saved_MAHOMES_II(features)

## use feature calculations to add any additional site metals to saved predictions 
save_preds = pd.read_csv("%s/features.csv"%(job_dir))
save_preds.set_index(['struc_id', 'metal1_resName', 'metal1_seqID'], inplace=True)

save_preds[['junk', 'site#']]=save_preds['SITE_ID'].str.rsplit("_", n=1, expand=True)
for col in ['metal2_resName', 'metal2_seqID', 'metal3_resName', 'metal3_seqID', 'metal4_resName', 'metal4_seqID']:
    if col not in save_preds.columns.tolist():
        save_preds[col]=""
    elif col in ['metal2_seqID', 'metal3_seqID', 'metal4_seqID']:
        save_preds[col]=save_preds[col].astype(int)

save_preds = save_preds[['site#', 'metal2_resName', 'metal2_seqID', 'metal3_resName', 'metal3_seqID', 'metal4_resName', 'metal4_seqID']].copy()
save_preds = save_preds.merge(preds[['prediction', 'bool_pred']].copy(), left_index=True, right_index=True)

## rename as needed
save_preds.reset_index(inplace=True)
save_preds.rename(columns={ 
     "struc_id":"jobid"
     ,"prediction":"percent catalytic predictions"
    }, inplace=True)
for i in range(1,5):
    save_preds.rename(columns={
     "metal%s_resName"%(str(i)):"Name%s"%(str(i))
     ,"metal%s_seqID"%(str(i)):"Res#%s"%(str(i))
    }, inplace=True)

save_preds['prediction']=""
save_preds.loc[save_preds["bool_pred"]==True,'prediction']="Catalytic"
save_preds.loc[save_preds["bool_pred"]==False,'prediction']="Not Catalytic"
save_preds.loc[save_preds["percent catalytic predictions"]=="",'prediction']="No prediction made"

save_preds["percent catalytic predictions"]=save_preds["percent catalytic predictions"]*100

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


save_preds.to_csv("%s/predictions.csv"%(job_dir), index=False)

