#!/bin/bash

#JOB_DIR="/home/mahomes/mahomes/test-MAHOMES-II-server"
JOB_DIR=$1
 
echo "Prepping structure files"
source ${MAHOMES_II_DIR}/venv/bin/activate
cd ${MAHOMES_II_DIR}
python prep_batch_job.py ${JOB_DIR}


## going over all input .pdbs and outputing the needed Rosetta, findgeo, bluues, and ghecom
echo "Starting third party calculations"
THIRD_PARTY_TOOLS=${MAHOMES_II_DIR}/FeatureCalculations/ThirdPartyCalculations
while read pdb; do
	input_file=${JOB_DIR}/${pdb}.pdb
	output_dir=${JOB_DIR}/${pdb}
	. ${THIRD_PARTY_TOOLS}/runParty3Calc.sh $pdb $input_file $output_dir
done < ${JOB_DIR}/batch_input.txt


## use the previous outputs to calculate the features for MAHOMES II input
echo "Working on calculating features"
cd ${MAHOMES_II_DIR}/FeatureCalculations
source ${MAHOMES_II_DIR}/venv/bin/activate
python batch_save_features.py ${JOB_DIR}

## make prediction(s)
cd ${MAHOMES_II_DIR}/MachineLearning/
source ${MAHOMES_II_DIR}/venv/bin/activate
python MakePredictions.py $JOB_DIR
cat $JOB_DIR/predictions.csv
