#!/bin/bash


JOB_DIR=$1

ROSETTA3="/var/www/apps/mahomes/rosetta/rosetta_src_2021.16.61629_bundle/main/source"
BLUUES_DIR="/var/www/apps/mahomes/bluues"
GHECOM_DIR="/var/www/apps/mahomes/ghecom/"
MAHOMES_II_DIR="/var/www/apps/mahomes/MAHOMES-II-server"

echo "Prepping structure files"
source ${MAHOMES_II_DIR}/venv/bin/activate
cd ${MAHOMES_II_DIR}
python prep_batch_job.py ${JOB_DIR}


## going over all input .pdbs and outputing the needed Rosetta, findgeo, bluues, and ghecom
echo "Starting third party calculations"
THIRD_PARTY_TOOLS=${MAHOMES_II_DIR}/FeatureCalculations/ThirdPartyCalculations
while read pdb; do
        start=`date +%s`
	input_file=${JOB_DIR}/${pdb}.pdb
	output_dir=${JOB_DIR}/${pdb}
	. ${THIRD_PARTY_TOOLS}/runParty3Calc.sh $pdb $input_file $output_dir $ROSETTA3 $BLUUES_DIR $GHECOM_DIR $MAHOMES_II_DIR
        end=`date +%s`
        echo "$pdb time: $((end-start))"
done < ${JOB_DIR}/batch_input.txt

start=`date +%s`
## use the previous outputs to calculate the features for MAHOMES II input
echo "Working on calculating features"
cd ${MAHOMES_II_DIR}/FeatureCalculations
source ${MAHOMES_II_DIR}/venv/bin/activate
python batch_save_features.py ${JOB_DIR}
end=`date +%s`
echo "Calculating features time: $((end-start))" 


start=`date +%s`
## make prediction(s)
cd ${MAHOMES_II_DIR}/MachineLearning/
source ${MAHOMES_II_DIR}/venv/bin/activate
python MakePredictions.py $JOB_DIR
cat $JOB_DIR/predictions.csv
end=`date +%s`
echo "predictions time: $((end-start))" 
