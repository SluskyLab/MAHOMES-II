import pandas as pd
# CV strats
from sklearn.model_selection import StratifiedKFold
#process results
from sklearn.metrics import confusion_matrix, matthews_corrcoef


# prep cv
def undersample(this_X, this_y, random_seed):
    y_Cat = this_y[this_y==True]
    y_nonCat = this_y[this_y==False]
    y_nonCat = y_nonCat#.sample(n=len(y_Cat)*3, axis=0, random_state=random_seed)
    y_return = pd.concat([y_Cat, y_nonCat])
    X_return = this_X.loc[y_return.index]
    return(X_return, y_return)

def check_aggregate_results(prediction_df):
    mcc = matthews_corrcoef(prediction_df['actual'], prediction_df['bool_pred'])
    TN, FP, FN, TP = confusion_matrix(prediction_df['actual'], prediction_df['bool_pred']).ravel()
    print(TP, TN, FP, FN)

    TPR=TP/(TP+FN)
    TNR=TN/(TN+FP)
    acc=(TP+TN)/(TP+TN+FP+FN)
    Prec=TP/(TP+FP)

    TPR*=100
    TNR*=100
    acc*=100
    Prec*=100
    #print(mcc, TPR, TNR, B_Acc, Acc, Prec)
    pub_df = pd.DataFrame([[acc, mcc, TPR, TNR, Prec,TP, TN, FP, FN]],
        columns=['Accuracy', 'MCC', 'Recall', 'TrueNegRate', 'Precision',"TP", "TN", "FP", "FN"])
    return(pub_df)


cv_type=StratifiedKFold(n_splits=7, random_state=None)
num_rand_seeds=10
def run_kfold(clf, sites_data, output=True):
    X = sites_data.loc[sites_data['Set']=='data'].copy()
    Y = X['Enzyme'].astype(bool); 
    del X['Enzyme']; del X['Set']
    print("X entries: %s \t features: %s"%(X.shape[0], X.shape[1]))

    # run CV
    kfold_results = []
    kfold_predictions=[]
    for i, (train_idx, test_idx) in enumerate(cv_type.split(X, Y)):
        X_other_folds, X_test_fold = X.iloc[train_idx].copy(), X.iloc[test_idx].copy()
        y_other_folds, y_test_fold = Y.iloc[train_idx].copy(), Y.iloc[test_idx].copy()
        if output:
            print("K-fold sizes (K=%s):"%(str(i)))
            print("\tfold: %s sites\t %s features"%(X_test_fold.shape[0], X_test_fold.shape[1]))
            print("\t\tnon-catalytic: %s \n\t\tcatalytic: %s"%(len(y_test_fold[y_test_fold==False]),len(y_test_fold[y_test_fold==True])))
            print("\ttraining: %s sites\t %s features"%(X_other_folds.shape[0], X_other_folds.shape[1]))
            print("\t\tnon-catalytic: %s \n\t\tcatalytic: %s"%(len(y_other_folds[y_other_folds==False]),len(y_other_folds[y_other_folds==True])))

        CV_site_preds = {'actual': pd.Series(y_test_fold, index=X_test_fold.index)}
        ## fit model to dataset and train on test set w/ diff random seeds
        for rand_seed in range(0,num_rand_seeds):
            X_train, y_train = undersample(X_other_folds, y_other_folds, rand_seed)
            if output:
                if rand_seed==0:
                    print("undersampling, random_state=%s, %s total sites"%(str(rand_seed), len(X_train)))
                    print("\t training set: \t%s non-catalytic sites \t %s catalytic sites"%(len(y_train[y_train==False]),len(y_train[y_train==True])))
            try:
                clf.set_params(random_state=rand_seed)
            except:
                print("no random state param")

            clf.fit(X_train.reset_index(drop=True), y_train.reset_index(drop=True))
            preds_CV = clf.predict(X_test_fold.reset_index(drop=True))
            CV_site_preds['prediction_%s'%(rand_seed)]= pd.Series(preds_CV, index=X_test_fold.index)

        CV_predictions = pd.DataFrame(CV_site_preds)
        CV_predictions['prediction']=0
        for rand_seed in range(0,num_rand_seeds):
            CV_predictions['prediction']+=CV_predictions['prediction_%s'%(rand_seed)]
        CV_predictions['prediction']=CV_predictions['prediction']/num_rand_seeds
        CV_predictions['bool_pred']=False
        CV_predictions.loc[CV_predictions['prediction']>=0.5, 'bool_pred']=True


        CV_scores = check_aggregate_results(CV_predictions)
        CV_scores['k']=i
        kfold_results.append(CV_scores)
        CV_predictions['k']=i
        kfold_predictions.append(CV_predictions)
    
    kfold_result_df = pd.concat(kfold_results, ignore_index=True)
    kfold_score = kfold_result_df[['Accuracy', 'MCC', 'Recall', 'TrueNegRate', 'Precision', 'TP', 'TN', 'FP', 'FN']].describe().loc[['mean']].reset_index(drop=True)
    kfold_std = kfold_result_df[['Accuracy', 'MCC', 'Recall', 'TrueNegRate', 'Precision', 'TP', 'TN', 'FP', 'FN']].describe().loc[['std']]
    for col in kfold_std.columns: kfold_score["%s_std"%(col)]=kfold_std[col].iloc[0]
    kfold_score = kfold_score[['Accuracy', 'Accuracy_std', 'MCC', 'MCC_std', 'Recall', 'Recall_std', 'Precision', 'Precision_std', 'TrueNegRate', 'TrueNegRate_std' ]].round(2)

    return(kfold_result_df, kfold_score, kfold_predictions)
        
def run_kfold_pipe(pipe, X, Y, output=True):
    print("X entries: %s \t features: %s"%(X.shape[0], X.shape[1]))

    # run CV
    kfold_results = []
    kfold_predictions=[]
    for i, (train_idx, test_idx) in enumerate(cv_type.split(X, Y)):
        X_other_folds, X_test_fold = X.iloc[train_idx].copy(), X.iloc[test_idx].copy()
        y_other_folds, y_test_fold = Y.iloc[train_idx].copy(), Y.iloc[test_idx].copy()
        if output:
            print("K-fold sizes (K=%s):"%(str(i)))
            print("\tfold: %s sites\t %s features"%(X_test_fold.shape[0], X_test_fold.shape[1]))
            print("\t\tnon-catalytic: %s \n\t\tcatalytic: %s"%(len(y_test_fold[y_test_fold==False]),len(y_test_fold[y_test_fold==True])))
            print("\ttraining: %s sites\t %s features"%(X_other_folds.shape[0], X_other_folds.shape[1]))
            print("\t\tnon-catalytic: %s \n\t\tcatalytic: %s"%(len(y_other_folds[y_other_folds==False]),len(y_other_folds[y_other_folds==True])))

        CV_site_preds = {'actual': pd.Series(y_test_fold, index=X_test_fold.index)}
        ## fit model to dataset and train on test set w/ diff random seeds
        for rand_seed in range(0,num_rand_seeds):
            X_train, y_train = undersample(X_other_folds, y_other_folds, rand_seed)
            if output:
                if rand_seed==0:
                    print("undersampling, random_state=%s, %s total sites"%(str(rand_seed), len(X_train)))
                    print("\t training set: \t%s non-catalytic sites \t %s catalytic sites"%(len(y_train[y_train==False]),len(y_train[y_train==True])))
            try:
                pipe.set_params(clf__random_state=rand_seed)
            except:
                print("no random state param")

            pipe.fit(X_train.reset_index(drop=True), y_train.reset_index(drop=True))
            preds_CV = pipe.predict(X_test_fold.reset_index(drop=True))
            CV_site_preds['prediction_%s'%(rand_seed)]= pd.Series(preds_CV, index=X_test_fold.index)

        CV_predictions = pd.DataFrame(CV_site_preds)
        CV_predictions['prediction']=0
        for rand_seed in range(0,num_rand_seeds):
            CV_predictions['prediction']+=CV_predictions['prediction_%s'%(rand_seed)]
        CV_predictions['prediction']=CV_predictions['prediction']/num_rand_seeds
        CV_predictions['bool_pred']=False
        CV_predictions.loc[CV_predictions['prediction']>=0.5, 'bool_pred']=True


        CV_scores = check_aggregate_results(CV_predictions)
        CV_scores['k']=i
        kfold_results.append(CV_scores)
        CV_predictions['k']=i
        kfold_predictions.append(CV_predictions)
    
    kfold_result_df = pd.concat(kfold_results, ignore_index=True)
    kfold_score = kfold_result_df[['Accuracy', 'MCC', 'Recall', 'TrueNegRate', 'Precision', 'TP', 'TN', 'FP', 'FN']].describe().loc[['mean']].reset_index(drop=True)
    kfold_std = kfold_result_df[['Accuracy', 'MCC', 'Recall', 'TrueNegRate', 'Precision', 'TP', 'TN', 'FP', 'FN']].describe().loc[['std']]
    for col in kfold_std.columns: kfold_score["%s_std"%(col)]=kfold_std[col].iloc[0]
    kfold_score = kfold_score[['Accuracy', 'Accuracy_std', 'MCC', 'MCC_std', 'Recall', 'Recall_std', 'Precision', 'Precision_std', 'TrueNegRate', 'TrueNegRate_std' ]].round(2)

    return(kfold_result_df, kfold_score, kfold_predictions)
 
