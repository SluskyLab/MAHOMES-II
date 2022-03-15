change
# MAHOMES II
Metal Activity Heuristic of Metalloprotein and Enzymatic Sites (MAHOMES) II - Predicts if a protein bound metal ion is enzymatic or non-enzymatic

## Overview
The ability to distinguish enzyme from non-enzyme sites remains an unsolved, challenging task. We've developed MAHOMES, a machine learning based tool which classifies metals bound to proteins as enzymatic or non-enzymatic. We intend to build on the previous work to make MAHOMES II, a more stable and robust version with a web server.

## System requirements
Feature calculations also require using Rosetta, FindGeo, and bluues which we run using Python 2.7 with CentOS 7 (or red hat now?).

## Installation guide
### set up virtual environment:
```
$ virtualenv --version # check for virtual enviroment
$ pip install virtualenv # download using pip (or pipx?) if no version is found
$ virtualenv -p /usr/bin/python3.8 venv # create new virtual environment
$ source venv/bin/activate # switch to new env (needs to be done whenever working on MAHOMES)
$ pip install -r requirements.txt # add packages to environment
```
Repeat above proccess for venv2 using python2.7 and requirements2.txt

## Instructions for use
This repository contains source code for ? different purposes

1. Making enzyme and non-enzyme predictions for a PDB with MAHOMES II
    After completing additional setup (FeatureCalculations/README.md and MachineLearning/README.md) the following steps can be used for making predicitons:
    1. add PDBs to ?
    2. run driver.sh
    3. The resulting data/<job_name>/sites_predictions.txt will contain the predictions for each site
        - final_prediction is MAHOMES enzyme or non-enzyme prediction
        - Prediction_<int> is the prediction made by one of the MAHOMES model using random seed <int>

2. Physicochemical feature calculations used by MAHOMES II (see FeatureCalculations/README.md for more information)

3. ML algorithm optimization
    MachineLearning/CV/MLwGrid.py can be ran to sue the nested cross validation to test the performance of an ML algorithm, feature set, and optmiization metric. Additional deatails and demo instructions are located at the top of the file.

4. MAHOMES training and performance evaluation
    This can be done by going threw each step of MachineLearning/MAHOMES_eval_T-metal-site.ipynb. This notebook includes:
    1. Reading in thee data-set and holdout T-metal-sites (followed by making the change to labels found during manual inspection)
    2. Scaling the features
    3. Under-sampling the training data and removing features that are not a part of the AllSph feature set
    4. Train an extra trees classifier to create a model that is makes predictions for T-metal-sites and is saved
    5. Repeat steps 3 and 4 with new random seeds nine times
    6. Averaging the predictions to get a final prediction which is rounded to 0 (non-enzyme) or 1 (enzyme)
    7. Calculating the performance metrics for the final T-metal-sites predictions
    This notebook should take less than five minutes.
    
    This can be demod using T-metal-sites for <job_name>. We have provided the calculated features and expected predictions. Note that saving new models may result in small differences due to the stochastic nature of random under-sampling the extra trees algorithm, but the final_prediction should remain the same. The demo should take less than one minute.


## Repo Contents

### Data

<b>publication_data</b> The data-set and T-metal-sites (sites.csv) and the calculated features for those sites

<b>T-metal-sites</b> T-metal-sites features and expected MAHOMES predictions.

### MachineLearning

The code used for the optimization, selection, and evaluation of the MAHOMES model

<b>CV</b>: Code used for nested cross validation which resulted in optimization and selection of MAHOMES

<b>MAHOMES_eval_T-metal-site.ipynb</b> Notebook that trains, evaluates, and saves the MAHOMES

<b>ModelsForMAHOMES</b> Saved versions of the scalers and extra trees models that make up MAHOMES

<b>MAHOMESNewPredictions.py</b> Uses saved version of MAHOMES to make predictions for sites with calculated features 
