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

##### import ########
# libraries
import numpy as np
import pandas as pd
import sys

##########################################
##########  data functions  ##############
##########################################
CAT_FEATS = ['geom_irr', 'geom_Reg', 'geom_Distort', 'geom_Filled', 'geom_PartFilled',
            'geom_cn2','geom_cn3','geom_cn4','geom_cn5','geom_cn6','geom_cn7','geom_cn8','geom_cn9']

def setDisplay(trainX, trainY, testX, testY):
    print("\nTRAIN entries: %s \t features: %s"%(trainX.shape[0], trainX.shape[1]))
    print("\tNum catalytic: %s \n\tNum non-catalytic: %s"%(len(trainY[trainY==1]),len(trainY[trainY==0])))
    print("CV entries: %s \t features: %s"%(testX.shape[0], testX.shape[1]))
    print("\tNum catalytic: %s \n\tNum non-catalytic: %s"%(len(testY[testY==1]),len(testY[testY==0])))

def get_data(feature_set, MLinput_dir = "../bin/ML_input", print_work=True):
    if 1 == feature_set: sites = pd.read_csv("%s/sites_feature_set1.csv"%(MLinput_dir))
    elif 2 == feature_set: sites = pd.read_csv("%s/sites_feature_set2.csv"%(MLinput_dir))
    elif 3 == feature_set: sites = pd.read_csv("%s/sites_feature_set3.csv"%(MLinput_dir))
    elif 4 == feature_set: sites = pd.read_csv("%s/sites_feature_set4.csv"%(MLinput_dir))
    elif 5 == feature_set: sites = pd.read_csv("%s/sites_feature_set5.csv"%(MLinput_dir))
    elif 6 == feature_set: sites = pd.read_csv("%s/sites_feature_set6.csv"%(MLinput_dir))
    else: 
        print("Invalid feature set")
        sys.exit()

    sites.set_index(['struc_id', 'metal1_resName', 'metal1_seqID'],inplace=True)

    ## get training/kfold sites, random under sample, and split out target value ("Catalytic")
    dataset_features = sites.loc[sites['Set']=='data'].copy()
    y_dataset = dataset_features["Enzyme"]
    X_dataset = dataset_features.copy()
    del X_dataset['Set']; del X_dataset['Enzyme']

    test_features = sites.loc[sites['Set']=='test'].copy()
    y_test = test_features["Enzyme"]
    X_test = test_features.copy()
    del X_test['Set']; del X_test['Enzyme']

    for term in CAT_FEATS:
        if term in X_dataset.columns:
            X_dataset[term]=X_dataset[term].astype('category')
            X_test[term]=X_test[term].astype('category')

    if print_work: setDisplay(X_dataset, y_dataset, X_test, y_test)
    return(X_dataset, y_dataset, X_test, y_test)
    
def get_data_corrected_labels(feature_set, MLinput_dir = "../bin/ML_input", print_work=True):
    if 1 == feature_set: sites = pd.read_csv("%s/sites_feature_set1.csv"%(MLinput_dir))
    elif 2 == feature_set: sites = pd.read_csv("%s/sites_feature_set2.csv"%(MLinput_dir))
    elif 3 == feature_set: sites = pd.read_csv("%s/sites_feature_set3.csv"%(MLinput_dir))
    elif 4 == feature_set: sites = pd.read_csv("%s/sites_feature_set4.csv"%(MLinput_dir))
    elif 5 == feature_set: sites = pd.read_csv("%s/sites_feature_set5.csv"%(MLinput_dir))
    elif 6 == feature_set: sites = pd.read_csv("%s/sites_feature_set6.csv"%(MLinput_dir))
    else: 
        print("Invalid feature set")
        sys.exit()

    sites.set_index(['struc_id', 'metal1_resName', 'metal1_seqID'],inplace=True)

    ## read in sites with corrected dataset and Tsite corrections found since more recenetly
    pipeline = pd.read_csv("../bin/MAHOMES_II_sites.csv", index_col=0)
    pipeline.set_index(['struc_id', 'resName1', 'seqNum1'],inplace=True)
    new_annotated_FPs = pd.read_excel("../bin/Tsite_FPs.xls", sheet_name="Tsite_FPs")

    ## remove recent Tsite pdb and sites according to manual annotations
    remove_pdbs =new_annotated_FPs.loc[new_annotated_FPs['annotation']=="remove", "pdb_name"].unique().tolist()
    remove_sites =new_annotated_FPs.loc[new_annotated_FPs['annotation']=="remove site", "meta_site"].unique().tolist()
    pipeline = pipeline.loc[~pipeline['pdb_name'].isin(remove_pdbs)]
    pipeline = pipeline.loc[~pipeline['meta_site'].isin(remove_sites)]
    
    ## the following failed feature calculations for MAHOMES or MAHOMES II, so remove to give same test-set
    remove_due_to_MHM_or_MHMII_failure = ['6h0f_C_0', '6jf1_A_1', '6izx_A_0', '6izx_A_1']
    pipeline = pipeline.loc[~pipeline['meta_site'].isin(remove_due_to_MHM_or_MHMII_failure)]
    sites=sites.loc[sites.index.isin(pipeline.index.tolist())].copy()

    ## fix recent Tsite label changes (and any in dataset that might have been different at time of feature calculation)
    change_labels =new_annotated_FPs.loc[new_annotated_FPs['annotation']=="Enzyme", "meta_site"].unique().tolist()
    pipeline.loc[pipeline['meta_site'].isin(change_labels), 'Enzyme']=True
    sites['Enzyme']=False
    sites.loc[sites.index.isin(pipeline.loc[pipeline['Enzyme']==True].index.tolist()),'Enzyme']=True

    ## get training/kfold sites, random under sample, and split out target value ("Catalytic")
    dataset_features = sites.loc[sites['Set']=='data'].copy()
    y_dataset = dataset_features["Enzyme"]
    X_dataset = dataset_features.copy()
    del X_dataset['Set']; del X_dataset['Enzyme']

    test_features = sites.loc[sites['Set']=='test'].copy()
    y_test = test_features["Enzyme"]
    X_test = test_features.copy()
    del X_test['Set']; del X_test['Enzyme']

    for term in CAT_FEATS:
        if term in X_dataset.columns:
            X_dataset[term]=X_dataset[term].astype('category')
            X_test[term]=X_test[term].astype('category')

    if print_work: setDisplay(X_dataset, y_dataset, X_test, y_test)
    return(X_dataset, y_dataset, X_test, y_test)
##########################################
##########  ML classifiers+  #############
##########################################
#classifiers
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression, RidgeClassifier, PassiveAggressiveClassifier
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, GradientBoostingClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis, LinearDiscriminantAnalysis

classifier_names = ["LogRegr", "Ridge", "PassAggr",
        "QDA", "LDA", "NaiveBayes",
        "NearNeigh",
        "LinSVM", "RBFSVM", "SigSVM",
        "RandomForest", "ExtraTrees", "GradBoost", 
        "NeurNet"]
parameter_space = [
    { "clf__penalty": ["l2", "l1"], "clf__C":[1.0, 0.01, 0.001], "clf__class_weight":['balanced', None]}, #LogRegr, removed l1
    { "clf__alpha": np.logspace(-4, 1, 6), "clf__class_weight":['balanced', None]}, #Ridge 
    { "clf__C": np.logspace(-1, 3, 5), "clf__class_weight":['balanced', None] }, #PassAggr
    { "clf__reg_param": np.linspace(0.5, 1, 7) }, #QDA
    { "clf__shrinkage": ["auto", 0, 0.1, 0.25, 0.5, 0.75, 1] }, #LDA
    { "clf__var_smoothing": np.logspace(-9, 0, 10) }, #Gauss
    [ {"clf__metric": ["minkowski"], "clf__p":[2, 3], "clf__n_neighbors": [5, 8, 10, 15]}, {"clf__metric": ["chebyshev"], "clf__n_neighbors": [5, 8, 10, 15]} ], #kNN
    { "clf__C": np.logspace(-3, 1, 5), "clf__class_weight":['balanced', None] }, #SVC lin
    { "clf__C": np.logspace(-3, 1, 5), "clf__gamma": ["scale", "auto"] }, #SVC rbf
    { "clf__C": np.logspace(-3, 1, 5),"clf__gamma": ["scale", "auto"] }, #SVC sig
    { "clf__class_weight": ["balanced", None], "clf__max_features": [0.25, 0.5, None], "clf__max_depth": [8, 16, 30] }, #RF
    { "clf__class_weight": ["balanced", None], "clf__max_features": [0.25, 0.5, None], "clf__max_depth": [8, 16, 30, 50, None] }, #ExtraTrees
    { "clf__learning_rate": [0.1, 0.005], "clf__max_features" : ['log2', 0.5, None], "clf__subsample":[0.5, 0.7, 0.9]},#GBClassifier
    { "clf__hidden_layer_sizes": [(50,), (100,), (200,)], "clf__alpha": [0.1, 0.01, 0.001] } #MLPClass
]
def get_classifier(classifier_name=None, classifier_index=None):
    ## idk why, but this takes for...ever if it is not within the function
    classifiers = [
    LogisticRegression(solver="liblinear"), 
    RidgeClassifier(),
    PassiveAggressiveClassifier(),
    QuadraticDiscriminantAnalysis(), 
    LinearDiscriminantAnalysis(solver = "lsqr"),
    GaussianNB(),
    KNeighborsClassifier(), 
    SVC(kernel="linear"), 
    SVC(kernel="rbf"), 
    SVC(kernel="sigmoid"),
    RandomForestClassifier(n_estimators=500, criterion='entropy', bootstrap=False),
    ExtraTreesClassifier(n_estimators=500, criterion='gini', bootstrap=False),
    GradientBoostingClassifier(criterion = 'friedman_mse', min_samples_split = 3, n_estimators=1000, loss='log_loss'),
    MLPClassifier(learning_rate_init = 0.01, activation='relu'),
    ]

    if (classifier_name!=None):
        classifier_index= classifier_names.index(classifier_name)
    elif (classifier_index!=None):
        classifier_name = classifier_names[classifier_index]
    else:
        print("failed to make valid classifier request")
        sys.exit()
    print(classifier_name, classifier_index)
    classifier = classifiers[classifier_index]
    clf_parameter_space = parameter_space[classifier_index]
    return(classifier_name, classifier, clf_parameter_space)

##########################################
##########      scaling     #############
##########################################
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler, RobustScaler, QuantileTransformer

scaler_names = ['Standard', 'Robust', 'QT_uniform', 'QT_normal',]
scalers = [StandardScaler(),
          RobustScaler(),
          QuantileTransformer(output_distribution="uniform"),
          QuantileTransformer(output_distribution="normal"),
     ]

def get_preprocessor(scale_features, catagorical_features, scaler_name=None, scaler_index=None):
    if (scaler_name!=None):
        scaler_index= scaler_names.index(scaler_name)
    elif (scaler_index!=None):
        scaler_name = scaler_names[scaler_index]
    else:
        print("failed to make valid scaler request")
        sys.exit()
    scales = scalers[scaler_index]
    numeric_preprocessor = Pipeline(steps=[
        ("impute_mean", SimpleImputer(strategy="mean")),
        ("scaler", scales),
        ])
    
    categorical_preprocessor = Pipeline(steps=[
            ("impute_median",SimpleImputer(strategy="median"),),
            #("minmax", MinMaxScaler()),
            ])
    preprocessor = ColumnTransformer([
            ("num", numeric_preprocessor, scale_features),
            ("cat", categorical_preprocessor, catagorical_features)],
            verbose_feature_names_out=False
            )
    return(scaler_name, preprocessor)
