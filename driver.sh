#!/bin/bash


THIRD_PARTY_TOOLS=${MAHOMES_II_DIR}/FeatureCalculations/ThirdPartyCalculations
JOB_DIR="/home/mahomes/mahomes/test-MAHOMES-II-server"

## I need to make the next three lines a loop, going over all input .pdbs as pdb_id
pdb_id="1"
input_file=${JOB_DIR}/${pdb_id}.pdb
output_dir=${JOB_DIR}/${pdb_id}
. ${THIRD_PARTY_TOOLS}/runParty3Calc.sh $pdb_id $input_file $output_dir


## make sure we are using the correct environment
echo "Working on calculating features"
source ${MAHOMES_II_DIR}/venv/bin/activate
cd ${MAHOMES_II_DIR}/FeatureCalculations
python batch_save_features.py ${JOB_DIR}

## make prediction(s)
cd ${MAHOMES_II_DIR}/MachineLearning/
python MakePredictions.py $JOB_DIR
cat $JOB_DIR/predictions.csv
