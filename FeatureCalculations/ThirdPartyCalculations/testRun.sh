#!/bin/bash

ROSETTA3="/path/to/rosetta/main/source"
BLUUES_DIR="/path/to/bluues"
GHECOM_DIR="/path/to/ghecom/"
MAHOMES_II_DIR="/path/to/MAHOMES_II"

THIRD_PARTY_TOOLS=${MAHOMES_II_DIR}/FeatureCalculations/ThirdPartyCalculations
pdb_id="1a7i_A"
echo $pdb_id
input_file=${THIRD_PARTY_TOOLS}/check_setup/${pdb_id}.pdb
output_dir=${THIRD_PARTY_TOOLS}/check_setup/output/
. ${THIRD_PARTY_TOOLS}/runParty3Calc.sh $pdb $input_file $output_dir $ROSETTA3 $BLUUES_DIR $GHECOM_DIR $MAHOMES_II_DIR

