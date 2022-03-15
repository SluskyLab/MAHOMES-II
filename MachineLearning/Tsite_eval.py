## This script trains an ML model using the metal ion site data-set and 
## evaluates the model's predictions on the T-metal-site set. 
## The models and scalers can be saved during this.

##### import ########
# libraries
import pandas as pd
import joblib
import os

#process results
from sklearn.metrics import confusion_matrix, matthews_corrcoef


## number of iterations to improve reproducability
num_rand_seeds = 10
def evaluate_model_with_Tsite(clf, data_df):
    X_train = data_df.loc[data_df['Set']=='data'].copy()
    Y_train = X_train['Enzyme'].astype(bool); 
    del X_train['Enzyme']; del X_train['Set']
    print("X train entries: %s \t features: %s"%(X_train.shape[0], X_train.shape[1]))

    X_test = data_df.loc[data_df['Set']=='test'].copy()
    Y_test = X_test['Enzyme'].astype(bool); 
    del X_test['Enzyme']; del X_test['Set']
    print("X test entries: %s \t features: %s"%(X_test.shape[0], X_test.shape[1]))

    feat_imp = []; raw_feat_imp = []
    ## get multiple predictions for test-set w/ diff random seeds
    test_site_preds = {'actual': pd.Series(Y_test, index=X_test.index)}
    for rand_seed in range(0,num_rand_seeds):
        ## train classifier and make test-set predictions
        clf.set_params(random_state=rand_seed)
        clf.fit(X_train, Y_train)
        test_preds = clf.predict(X_test)
        test_site_preds['prediction_%s'%(rand_seed)]= pd.Series(test_preds, index=X_test.index)
        
        ## output results for this random seed to get an idea of prediction variation levels
        TN, FP, FN, TP = confusion_matrix(Y_test, test_preds).ravel()
        mcc = matthews_corrcoef(Y_test, test_preds)
        print("\tTP=%s \tTN=%s \tFP=%s \tFN=%s"%(TP, TN, FP, FN))
        
        new_feat_imp = pd.DataFrame([clf.feature_importances_], columns=X_train.columns)
        max_feat_imp =new_feat_imp.iloc[0].max()
        new_feat_imp.iloc[0] = (new_feat_imp.iloc[0]/max_feat_imp)*100
        feat_imp.append(new_feat_imp)
        raw_feat_imp.append(pd.DataFrame([clf.feature_importances_], columns=X_train.columns))


    ## calcualte the average of all random seed predictions
    test_predictions = pd.DataFrame(test_site_preds)
    test_predictions['prediction']=0
    for rand_seed in range(0,num_rand_seeds):
        test_predictions['prediction']+=test_predictions['prediction_%s'%(rand_seed)] 
    test_predictions['prediction']=test_predictions['prediction']/num_rand_seeds
    
    ## make final prediction
    test_predictions['bool_pred']=False
    test_predictions.loc[test_predictions['prediction']>=0.5, 'bool_pred']=True
    raw_feat_imp_df = pd.concat(raw_feat_imp, ignore_index=True)
    feat_imp_df = pd.concat(feat_imp, ignore_index=True)
    total_imp_df = [];
    for term in feat_imp_df.columns:
        total_imp = feat_imp_df[term].sum()/num_rand_seeds
        total_imp_df.append([term, total_imp])
    total_imp_df = pd.DataFrame(total_imp_df, columns=['Feature', 'Importance'])
    total_imp_df=total_imp_df.sort_values('Importance', ascending=False).reset_index(drop=True)

    return(test_predictions, total_imp_df, raw_feat_imp_df)



def evaluate_pipe_with_Tsite(pipe, X_train, Y_train, X_test, Y_test, return_feature_imp=False, print_work=True, save_models=False, save_dir="saved_models"):
    if print_work:
        print("running evaluate_pipe_with_Tsite")
        print("X train entries: %s \t features: %s"%(X_train.shape[0], X_train.shape[1]))
        print("X test entries: %s \t features: %s"%(X_test.shape[0], X_test.shape[1]))
    feat_imp = []; raw_feat_imp = []
    ## get multiple predictions for test-set w/ diff random seeds
    test_site_preds = {'actual': pd.Series(Y_test, index=X_test.index)}
    for rand_seed in range(0,num_rand_seeds):
        ## train classifier and make test-set predictions
        try:
            pipe.set_params(clf__random_state=rand_seed)
        except:
            if print_work: print("just repeating w/o changing anything :)")
        pipe.fit(X_train, Y_train)
        if save_models==True:
            os.system("mkdir -p %s"%save_dir)
            joblib.dump(pipe, "%s/MAHOMES_II_r%s.pkl"%(save_dir, rand_seed))

        test_preds = pipe.predict(X_test)
        test_site_preds['prediction_%s'%(rand_seed)]= pd.Series(test_preds, index=X_test.index)
        
        ## output results for this random seed to get an idea of prediction variation levels
        TN, FP, FN, TP = confusion_matrix(Y_test, test_preds).ravel()
        if print_work: print("\t%s: TP=%s \tTN=%s \tFP=%s \tFN=%s"%(rand_seed, TP, TN, FP, FN))
        if return_feature_imp:
            try:
                new_feat_imp = pd.DataFrame([pipe.steps[1][1].feature_importances_], columns=X_train.columns)
                max_feat_imp =new_feat_imp.iloc[0].max()
                new_feat_imp.iloc[0] = (new_feat_imp.iloc[0]/max_feat_imp)*100
                feat_imp.append(new_feat_imp)
                raw_feat_imp.append(pd.DataFrame([pipe.steps[1][1].feature_importances_], columns=X_train.columns))
            except:
                print("could not get model feature imp.")

    ## calcualte the average of all random seed predictions
    test_predictions = pd.DataFrame(test_site_preds)
    test_predictions['prediction']=0
    for rand_seed in range(0,num_rand_seeds):
        test_predictions['prediction']+=test_predictions['prediction_%s'%(rand_seed)] 
    test_predictions['prediction']=test_predictions['prediction']/num_rand_seeds
    ## make final prediction
    test_predictions['bool_pred']=False
    test_predictions.loc[test_predictions['prediction']>=0.5, 'bool_pred']=True
    
    struc_preds = flatten_to_structure_predictions(test_predictions, pipe_file = "../bin/MAHOMES_II_sites.csv", print_work=print_work)
    struc_score = check_result_metrics(struc_preds, print_work=print_work)
    struc_score= structure_pred_reproducability(struc_preds, struc_score)


    if return_feature_imp and (len(raw_feat_imp)>0):
        raw_feat_imp_df = pd.concat(raw_feat_imp, ignore_index=True)
        feat_imp_df = pd.concat(feat_imp, ignore_index=True)
        total_imp_df = [];
        for term in feat_imp_df.columns:
            total_imp = feat_imp_df[term].sum()/num_rand_seeds
            total_imp_df.append([term, total_imp])
        total_imp_df = pd.DataFrame(total_imp_df, columns=['Feature', 'Importance'])
        total_imp_df=total_imp_df.sort_values('Importance', ascending=False).reset_index(drop=True)
    else:
        total_imp_df=[]; raw_feat_imp_df=[]

    return(test_predictions, struc_preds, struc_score, total_imp_df, raw_feat_imp_df)


## return result metrics for final predictions
def check_result_metrics(prediction_df, print_work=True, include_convergence_score=True):
    mcc = matthews_corrcoef(prediction_df['actual'], prediction_df['bool_pred'])
    TN, FP, FN, TP = confusion_matrix(prediction_df['actual'], prediction_df['bool_pred']).ravel()
    if print_work: print("TN %s, FP %s, FN %s, TP %s"%(TN, FP, FN, TP))
    TPR=(TP/(TP+FN))*100
    TNR=(TN/(TN+FP))*100
    acc=((TP+TN)/(TP+TN+FP+FN))*100
    Prec=(TP/(TP+FP))*100
    if include_convergence_score:
        convergence_score =(prediction_df['convergence'].sum()/len(prediction_df))*100
        return(pd.DataFrame([[acc, mcc, TPR, Prec, TNR, convergence_score]],
            columns=['Accuracy', 'MCC', 'Recall', 'Precision', 'TrueNegRate', 'convergence score']))
    return(pd.DataFrame([[acc, mcc, TPR, Prec, TNR]],
            columns=['Accuracy', 'MCC', 'Recall', 'Precision', 'TrueNegRate']))
    
def flatten_to_structure_predictions(Tsite_predictions, pipe_file = "../bin/MAHOMES_II_sites.csv", print_work=True):
    Tsite_predictions=Tsite_predictions.reset_index()
    testset_sites = pd.read_csv(pipe_file, index_col=0)
    testset_sites = testset_sites.loc[testset_sites['Set']=='test']
    testset = pd.merge(Tsite_predictions, testset_sites, how='inner', left_on=['struc_id', 'metal1_resName', 'metal1_seqID'], right_on=['struc_id', 'resName1', 'seqNum1'])

    struc_preds = []
    for meta_site_id in testset['meta_site'].unique():
        meta_site_preds = testset.loc[testset['meta_site']==meta_site_id]

        ## calculate this sites convergence
        confergence=0
        for i, row in meta_site_preds.iterrows():
            if row['bool_pred']:confergence+=1
            else: confergence-=1
        confergence = abs(confergence/len(meta_site_preds))

        ## calculate this sites final predictions
        if len(meta_site_preds)%10!=0:
            if print_work: print("not 10 repeat meta-site:", meta_site_id)
        total_pred = meta_site_preds['bool_pred'].sum()
        num_pred = len(meta_site_preds)
        pred = total_pred/num_pred
        bool_pred = False
        if pred>=0.5:
            bool_pred=True
        if len(meta_site_preds['actual'].unique())<1:
            print(meta_site_id)
        actual = meta_site_preds['actual'].unique().tolist()[0]
        if len(meta_site_preds['actual'].unique())>1:
            print("F!", meta_site_id, meta_site_preds['actual'].unique().tolist())
            #display(meta_site_preds)
        this_pred_df = pd.DataFrame([[meta_site_id, total_pred, num_pred, pred, bool_pred, actual, confergence]]
                                    , columns=['meta_site_id', 'total_pred', 'num_pred', 'pred', 'bool_pred', 'actual', 'convergence'])
        struc_preds.append(this_pred_df)
    struc_preds = pd.concat(struc_preds, ignore_index=True)
    return(struc_preds)


def structure_pred_reproducability(struc_preds, score_df):
    #print("%s predictions"%len(struc_preds))
    not_enzyme_preds = struc_preds['pred']==0
    score_df['not_enzyme_preds'] = len(struc_preds.loc[not_enzyme_preds])
                           
    enzyme_preds = struc_preds['pred']==1
    score_df['enzyme_preds'] = len(struc_preds.loc[enzyme_preds])
                           
    score_df['other_preds'] =len( struc_preds.loc[~enzyme_preds & ~not_enzyme_preds])

    bad_pred1 = struc_preds['pred']<=0.75
    bad_pred2 = struc_preds['pred']>=0.25
    score_df['bad_preds'] = len(struc_preds.loc[bad_pred1 & bad_pred2])
    return(score_df)


def output_structure_pred_reproducability(struc_preds):
    print("%s predictions"%len(struc_preds))
    not_enzyme_preds = struc_preds['pred']==0
    print("%s complete non-enzyme predictions"%len(struc_preds.loc[not_enzyme_preds]))
    enzyme_preds = struc_preds['pred']==1
    print("%s complete enzyme predictions"%len(struc_preds.loc[enzyme_preds]))
    print("%s non-complete predictions"%len(struc_preds.loc[~enzyme_preds & ~not_enzyme_preds]))

    bad_pred1 = struc_preds['pred']<=0.75
    bad_pred2 = struc_preds['pred']>=0.25
    print("%s predictions from 0.25 to 0.75"%len(struc_preds.loc[bad_pred1 & bad_pred2]))




