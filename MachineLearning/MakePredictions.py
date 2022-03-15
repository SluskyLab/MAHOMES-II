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
for col in ['metal2_resName', 'metal2_seqID', 'metal3_resName', 'metal3_seqID', 'metal4_resName', 'metal4_seqID']:
    if col not in save_preds.columns.tolist():
        save_preds[col]=""

save_preds = save_preds[['metal2_resName', 'metal2_seqID', 'metal3_resName', 'metal3_seqID', 'metal4_resName', 'metal4_seqID']].copy()
save_preds = save_preds.merge(preds[['prediction', 'bool_pred']].copy(), left_index=True, right_index=True)

## rename as needed
save_preds.reset_index(inplace=True)
save_preds.rename(columns={ 
     "struc_id":"jobid"
    ,"bool_pred":"enzyme"
    }, inplace=True)
for i in range(1,5):
    save_preds.rename(columns={
     "metal%s_resName"%(str(i)):"resName%s"%(str(i))
     ,"metal%s_seqID"%(str(i)):"resNum%s"%(str(i))
    }, inplace=True)

## add user input file name
with open("%s/input-legend.json"%(job_dir), "r") as f:
    name_dict = json.loads(f.read())
save_preds['jobid']=save_preds['jobid']+".pdb"
save_preds['user_input']=save_preds['jobid'].map(name_dict)

## reorder, which will change when I change col names!
save_preds = save_preds[[
    'jobid','user_input'
    ,'resName1', 'resNum1', 'resName2', 'resNum2' #,'metal1_resName', 'metal1_seqID', 'metal2_resName', 'metal2_seqID'
    ,'resName3', 'resNum3', 'resName4', 'resNum4' #,'metal3_resName', 'metal3_seqID', 'metal4_resName', 'metal4_seqID'
    ,'prediction', 'enzyme']].copy()

save_preds.to_csv("%s/predictions.csv"%(job_dir))

