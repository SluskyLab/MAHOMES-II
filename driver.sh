#! /bin/bash
JOB_NAME=$1
JOB_DIR="/var/www/apps/mahomes/test_data_for_mahomes2021-server/${JOB_NAME}"


# use python3.8 on server from venv,
# if that fails, work on adding MLenv from cluster to server
source venv/bin/activate
cd FeatureCalculations/
FEAT_DIR=$(pwd)
#echo $FEAT_DIR

## identify all sites on input structures and create input files for batch feature calculations
python PrepSITES.py $JOB_DIR
echo "Prep complete"

## create *_ByResRelaxValues.txt from input PDB file
cd $FEAT_DIR
./BatchScore.sh $JOB_DIR
echo "Scoring complete"

## maybe TODO: fix grid/vol naming. It is currently the only way that will run,
##         so figure out what needs to be changed and the ideal thing to change it to
## Use Rosetta pocket_grid to generate info for pocket and lining
cd $FEAT_DIR
./BatchGrid.sh $JOB_DIR
echo "pocket_grid complete"

## maybe TODO: pdb2pqr has a lot of warnings, which I think are mostly related to reading in the Rosetta output pdb, but look into this
## use pdb2pqr and bluues to generate *_ElectFeatures
cd $FEAT_DIR
./BatchTitr.sh $JOB_DIR
echo "Bluues complete"

## make  *.findgeo output files
cd $FEAT_DIR
./BatchGeom.sh $JOB_DIR
echo "Geom complete"

## combine all features into one file
cd $FEAT_DIR
source ../venv/bin/activate
python WriteFeaturesToList.py $JOB_DIR


cd $FEAT_DIR
cd ../
echo "finished features section of driver"

## make prediction(s)
cd MachineLearning/
python MAHOMESNewPredictions.py $JOB_DIR
cd ..
cat $JOB_DIR/sites_predictions.txt

