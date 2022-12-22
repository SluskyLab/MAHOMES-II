import pandas as pd
import Tsite_eval as Tsite
import MAHOMES_II as MHMII
import warnings
warnings.filterwarnings("ignore")


MHMII_alg = "GradBoost"
MHMII_scl = "QT_uniform"
MHMII_FeatureSet=4
MHMII_params="{'clf__learning_rate': 0.05, 'clf__max_features': 'log2', 'clf__subsample': 0.50, 'clf__n_estimators': 5000, 'clf__max_depth': 6}"
MHMII.run_and_save_MAHOMES_II(MHMII_alg, MHMII_scl, MHMII_FeatureSet, MHMII_params)


AFset_features = pd.read_csv("../bin/AFset_features.csv")
preds = MHMII.make_predictions_with_saved_MAHOMES_II(AFset_features)
print(preds.shape)

AFset_preds = preds.copy()
print(AFset_preds.shape)
AFset_sites = pd.read_csv("../bin/AF_sites.csv", index_col=0)
AFset_sites.set_index(['struc_id', 'metal1_resName', 'metal1_seqID'], inplace=True)
AFset_preds = AFset_preds.merge(AFset_sites, left_index=True, right_index=True)
AFset_preds = AFset_preds.loc[AFset_preds['prediction'].notnull()]
AFset_preds['actual']=AFset_preds['Enzyme']
AFset_preds['min residue pLDDT within 5Å of site']=AFset_preds['min_pLDDT_5'].astype(float)
AFset_preds['min residue pLDDT within 10Å of site']=AFset_preds['min_pLDDT_10'].astype(float)
AFset_preds['min residue pLDDT within 15Å of site']=AFset_preds['min_pLDDT_15'].astype(float)
print(AFset_preds.shape)

all_metrics = []
metrics = Tsite.check_result_metrics(AFset_preds, print_work=True, include_divergence_metrics=False)
sizes = AFset_preds.groupby(['Enzyme']).size()
metrics['enzyme sites']= sizes[1]; metrics['not-enzyme sites']= sizes[0]; metrics['total sites']= AFset_preds.shape[0]
metrics['Set']="All made AF predictions"
all_metrics.append(metrics)

for col in ['min residue pLDDT within 5Å of site', 'min residue pLDDT within 10Å of site', 'min residue pLDDT within 15Å of site']:
    AFset_preds = AFset_preds.loc[AFset_preds[col]>=90]
    metrics = Tsite.check_result_metrics(AFset_preds, print_work=True, include_convergence_score=False)
    sizes = AFset_preds.groupby(['Enzyme']).size()
    metrics['enzyme sites']= sizes[1]; metrics['not-enzyme sites']= sizes[0]; metrics['total sites']= AFset_preds.shape[0];
    metrics['Set']="%s > 90"%col
    all_metrics.append(metrics)

all_metrics=pd.concat(all_metrics, ignore_index=True)
print(all_metrics[['Set', 'enzyme sites', 'not-enzyme sites', 'Accuracy', 'MCC', 'Recall', 'TrueNegRate', 'Precision']].round(2))
