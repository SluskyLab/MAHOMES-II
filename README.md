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
$ pip install virtualenv # download using pip if no version is found
$ virtualenv -p /usr/bin/python3 venv # create new virtual environment
$ source venv/bin/activate # switch to new env 
$ pip install -r requirements.txt # add packages to environment
```
Repeat above proccess for venv2 using python2.7 and requirements2.txt

## Follow FeatureCalculations/README.md to set up third-party feature calculations

## Training and saving MAHOMES II ML models
```
$ ./prep_ML.sh
```

## Making enzyme and non-enzyme predictions for a PDB with MAHOMES II
1. Make a new directory and place one or more PDB files in it
2. Run the following command, replacing $JOB_DIR with the directory from step 1
```
. /path-to-repo/MAHOMES-II-server/driver.sh $JOB_DIR
```
The directory should now have predictions.csv as well as all the calculated features.

