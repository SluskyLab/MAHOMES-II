import pandas as pd
import GeneralML as ML
import kfold_eval as kfold
import Tsite_eval as Tsite
from sklearn.pipeline import make_pipeline, Pipeline
import ast
import joblib

MAHOMES_II_feature_order = [
    ## electrostatics inner shell
    'Elec_ins_Z_mu2_count', 'Elec_ins_max_DestabRank', 'Elec_ins_max_GBR6', 'Elec_ins_max_ResSigDev', 'Elec_ins_max_SolvEnergy', 'Elec_ins_max_SolvExp', 'Elec_ins_max_StabRank', 'Elec_ins_max_dpKa_bg', 'Elec_ins_max_dpKa_desolv', 'Elec_ins_max_dpKa_titr', 'Elec_ins_max_mu2', 'Elec_ins_max_mu4', 'Elec_ins_max_pKa_shift'
    ,'Elec_ins_mean_dpKa_desolv', 'Elec_ins_mean_dpKa_titr', 'Elec_ins_mean_pKa_shift', 'Elec_ins_min_DestabRank', 'Elec_ins_min_ResSigDev', 'Elec_ins_min_SolvEnergy', 'Elec_ins_min_SolvExp', 'Elec_ins_min_StabRank', 'Elec_ins_min_dpKa_bg', 'Elec_ins_min_dpKa_desolv', 'Elec_ins_min_dpKa_titr', 'Elec_ins_min_mu2'
    ,'Elec_ins_range_SolvEnergy', 'Elec_ins_range_SolvExp', 'Elec_ins_range_dpKa_bg', 'Elec_ins_range_dpKa_titr', 'Elec_ins_range_mu2', 'Elec_ins_range_pKa_shift'
    ## electrostatics outer shell
    ,'Elec_outs_Z_mu2_count', 'Elec_outs_max_DestabRank', 'Elec_outs_max_GBR6', 'Elec_outs_max_ResSigDev', 'Elec_outs_max_SolvEnergy', 'Elec_outs_max_SolvExp', 'Elec_outs_max_StabRank', 'Elec_outs_max_dpKa_bg', 'Elec_outs_max_dpKa_titr', 'Elec_outs_max_mu2', 'Elec_outs_max_pKa_shift'
    ,'Elec_outs_mean_DestabRank', 'Elec_outs_mean_SolvEnergy', 'Elec_outs_mean_SolvExp', 'Elec_outs_mean_StabRank', 'Elec_outs_mean_dpKa_bg', 'Elec_outs_mean_dpKa_desolv', 'Elec_outs_mean_dpKa_titr', 'Elec_outs_mean_mu2', 'Elec_outs_mean_pKa_shift'
    ,'Elec_outs_min_DestabRank', 'Elec_outs_min_GBR6', 'Elec_outs_min_ResSigDev', 'Elec_outs_min_SolvEnergy', 'Elec_outs_min_SolvExp', 'Elec_outs_min_StabRank', 'Elec_outs_min_dpKa_bg', 'Elec_outs_min_dpKa_desolv', 'Elec_outs_min_dpKa_titr', 'Elec_outs_min_mu3', 'Elec_outs_min_pKa_shift'
    ,'Elec_outs_numResis', 'Elec_outs_range_SolvExp', 'Elec_outs_range_dpKa_bg', 'Elec_outs_range_mu3', 'Elec_outs_std_dpKa_titr', 'Elec_outs_std_pKa_shift'
    ## electrostatics xenv
    ,'Elec_xenv_DestabRank', 'Elec_xenv_SolvExp', 'Elec_xenv_StabRank', 'Elec_xenv_mu4'
    ## geom
    ,'MetalAtoms', 'SiteAtoms', 'geom_AtomRMSD', 'geom_AvgN', 'geom_AvgO', 'geom_AvgOther', 'geom_AvgS', 'geom_Distort', 'geom_Filled', 'geom_LigO', 'geom_PartFilled', 'geom_Reg', 'geom_cn2', 'geom_cn3', 'geom_cn4', 'geom_cn5', 'geom_cn6', 'geom_cn7', 'geom_irr'
    ## New pocket features
    ,'pkt_info_ClusterNum', 'pkt_info_Ngrid', 'pkt_info_Rinac(A) av', 'pkt_info_Rinac(A) mi', 'pkt_info_SITE_pocket_distance_min', 'pkt_info_metal_depth', 'pkt_info_metal_height', 'pkt_info_num_adjacent_pockets', 'pkt_info_pocket_depth'
    ## Pocket lining
    ,'pkt_lining_SC_vol_perc', 'pkt_lining_avg_eisen_hp', 'pkt_lining_avg_kyte_hp', 'pkt_lining_max_eisen', 'pkt_lining_min_eisen', 'pkt_lining_min_kyte', 'pkt_lining_skew_eisen', 'pkt_lining_skew_kyte', 'pkt_lining_std_dev_eisen', 'pkt_lining_std_dev_kyte'
    ## Pocket void
    ,'pkt_shape_Depth', 'pkt_shape_LongAxHigh', 'pkt_shape_LongAxLow', 'pkt_shape_LongPath', 'pkt_shape_OffsetHigh', 'pkt_shape_OffsetLow', 'pkt_shape_OffsetMid', 'pkt_shape_SITEDistCenter', 'pkt_shape_SITEDistNormCenter', 'pkt_shape_ShortAxMid'
    ## Rosetta energy terms
    ,'rbf_dslf_fa13', 'rbf_fa_atr', 'rbf_fa_dun', 'rbf_fa_dun_dev', 'rbf_fa_dun_rot', 'rbf_fa_elec', 'rbf_fa_intra_atr_xover4', 'rbf_fa_intra_elec', 'rbf_fa_intra_rep_xover4', 'rbf_fa_intra_sol_xover4', 'rbf_fa_rep'
    ,'rbf_fa_sol', 'rbf_hbond_bb_sc', 'rbf_hbond_sc', 'rbf_hxl_tors', 'rbf_lk_ball_bridge', 'rbf_numScoredResis', 'rbf_omega', 'rbf_p_aa_pp', 'rbf_pro_close', 'rbf_rama_prepro', 'rbf_ref', 'rbf_total', 'rbf_yhh_planarity'
    ]


def set_pipe_clf_params(params, pipe, print_params=False):
     use_params = ast.literal_eval(params)
     pipe = pipe.set_params(**use_params)
     if print_params:
          cur_params=pipe.named_steps['clf'].get_params()
          for param in cur_params:
               print(param, cur_params[param])
     return(pipe)

def eval_kfold(X_train, Y_train, clf_name, preprocess, use_params):
    _, clf, _ = ML.get_classifier(classifier_name=clf_name)
    pipe = Pipeline(steps=[("preprocess", preprocess), ("clf", clf)])
    
    this_pipe = set_pipe_clf_params(use_params, pipe)
    kfold_result_df, kfold_score, kfold_predictions = kfold.run_kfold_pipe(this_pipe, X_train, Y_train, output=False)
    print("CV:")
    print(kfold_score.round(2))
    return(kfold_result_df, kfold_score, kfold_predictions)

def eval_Tsite(X_train, Y_train, X_test, Y_test, clf_name, preprocess, use_params, save_mdls_dir):
    _, clf, _ = ML.get_classifier(classifier_name=clf_name)
    pipe = Pipeline(steps=[("preprocess", preprocess), ("clf", clf)])
    this_pipe = set_pipe_clf_params(use_params, pipe, print_params=False)

    all_predictions, struc_preds, struc_score, _, _ = Tsite.evaluate_pipe_with_Tsite(this_pipe, X_train, Y_train, X_test, Y_test, print_work=True, save_models=True, save_dir=save_mdls_dir)
    
    return_score=struc_score.copy()
    for col in struc_score.columns:
         return_score["Tsites_%s"%col]=return_score[col]; del return_score[col]
    print("Tsites:")
    print(return_score.round(2))
    return(return_score)

def run_and_save_MAHOMES_II(this_classifier_name, this_scaler_name, this_feature_set, this_use_params):
        
    X_train, Y_train, X_test, Y_test  = ML.get_data_corrected_labels(feature_set=this_feature_set, print_work=True)
    X_train=X_train[MAHOMES_II_feature_order]; X_test=X_test[MAHOMES_II_feature_order]
    print(X_train.shape, X_test.shape)
    
    clf_name, clf, prm_space = ML.get_classifier(classifier_name=this_classifier_name)
    cont_features =X_train.select_dtypes(exclude=['category']).columns
    print("num cont feats = %s"%(len(cont_features)))
    catg_features =X_train.select_dtypes(['category']).columns
    print("num catg feats = %s"%(len(catg_features)))
    preprocess_name, preprocess = ML.get_preprocessor(cont_features, catg_features, scaler_name=this_scaler_name)
    eval_df = pd.DataFrame(["MAHOMES II"], columns=['model'])
    
    this_kfold_result_df, this_kfold_score, this_kfold_predictions= eval_kfold(X_train, Y_train, clf_name, preprocess, this_use_params)
    for stat in ['TP', 'TN', 'FP', 'FN']:
        print("%s:%s"%(stat, this_kfold_result_df[stat].sum()), end='\t\t')
        eval_df[stat]=this_kfold_result_df[stat].sum()
    print('\n')
    for stat in ['MCC', 'Precision', 'Recall', 'TrueNegRate', 'Accuracy']:
        print("%s\t%s\t%s"%(stat, round(this_kfold_result_df[stat].mean(), 4), round(this_kfold_result_df[stat].std(), 4)))
        eval_df[stat]=this_kfold_result_df[stat].mean()
        eval_df["%s_std"%stat]=this_kfold_result_df[stat].std()

    save_dir="saved_models"
    tsites_return_score = eval_Tsite(X_train, Y_train, X_test, Y_test, clf_name, preprocess, this_use_params, save_mdls_dir=save_dir)

    for term in tsites_return_score.columns:
        eval_df[term]=tsites_return_score[term]
    eval_df.to_csv("%s/MAHOMES_II_build_results.csv"%(save_dir))
    print(eval_df.round(2))

#MHMII_alg = "GradBoost"
#MHMII_scl = "QT_uniform"
#MHMII_FeatureSet=4
#MHMII_params="{'clf__learning_rate': 0.05, 'clf__max_features': 'log2', 'clf__subsample': 0.50, 'clf__n_estimators': 5000, 'clf__max_depth': 6}"
#run_and_save_MAHOMES_II(MHMII_alg, MHMII_scl, MHMII_FeatureSet, MHMII_params)

num_rand_seeds=10
def make_predictions_with_saved_MAHOMES_II(X_new, save_dir="saved_models"):
    X_new.set_index(['struc_id', 'metal1_resName', 'metal1_seqID', 'metal1_serial'],inplace=True)
    new_site_preds = {'prediction': pd.Series("-1", index=X_new.index)}

    X_dont_predict =X_new.loc[X_new['bad_site']==True].copy()
    X_new=X_new.loc[X_new['bad_site']==False].copy()
    
    for term in ML.CAT_FEATS:
        if term in X_new.columns:
            X_new[term]=X_new[term].astype('category')
    X_new=X_new[MAHOMES_II_feature_order]

    ## get multiple predictions for test-set w/ diff random seeds
    #new_site_preds = {'prediction': pd.Series("-1", index=X_new.index)}
    for rand_seed in range(0,num_rand_seeds):
        pipe= joblib.load( "%s/MAHOMES_II_r%s.pkl"%(save_dir, rand_seed))
        new_preds = pipe.predict(X_new)
        new_site_preds['prediction_%s'%(rand_seed)]= pd.Series(new_preds, index=X_new.index)

    ## calcualte the average of all random seed predictions
    new_predictions = pd.DataFrame(new_site_preds)
    new_predictions['prediction']=0
    for rand_seed in range(0,num_rand_seeds):
        new_predictions['prediction']+=new_predictions['prediction_%s'%(rand_seed)] 
    new_predictions['prediction']=new_predictions['prediction']/num_rand_seeds
    ## make final prediction
    new_predictions['bool_pred']=False
    new_predictions.loc[new_predictions['prediction']>=0.5, 'bool_pred']=True
    return(new_predictions)


