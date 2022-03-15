#!/bin/bash


THIRD_PARTY_TOOLS=${MAHOMES_II_DIR}/FeatureCalculations/ThirdPartyCalculations
#JOB_DIR="/home/mahomes/mahomes/test-MAHOMES-II-server"
JOB_DIR=$1


## 
echo "Prepping structure files"
source ${MAHOMES_II_DIR}/venv/bin/activate
cd ${MAHOMES_II_DIR}
python prep_batch_job.py ${JOB_DIR}


## I need to make the next three lines a loop, going over all input .pdbs as pdb_id
echo "Starting third party calculations"
while read pdb; do
	input_file=${JOB_DIR}/${pdb}.pdb
	output_dir=${JOB_DIR}/${pdb}
	#echo "$pdb"
	#echo "$input_file"
	#echo "$output_dir"
	. ${THIRD_PARTY_TOOLS}/runParty3Calc.sh $pdb $input_file $output_dir
done < ${JOB_DIR}/batch_input.txt

#pdb_id="1"
#input_file=${JOB_DIR}/${pdb_id}.pdb
#output_dir=${JOB_DIR}/${pdb_id}
#. ${THIRD_PARTY_TOOLS}/runParty3Calc.sh $pdb_id $input_file $output_dir


## make sure we are using the correct environment
echo "Working on calculating features"
source ${MAHOMES_II_DIR}/venv/bin/activate
cd ${MAHOMES_II_DIR}/FeatureCalculations
python batch_save_features.py ${JOB_DIR}

## make prediction(s)
cd ${MAHOMES_II_DIR}/MachineLearning/
python MakePredictions.py $JOB_DIR
cat $JOB_DIR/predictions.csv
