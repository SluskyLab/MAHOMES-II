#!/bin/bash

THIRD_PARTY_TOOLS=${MAHOMES_II_DIR}/FeatureCalculations/ThirdPartyCalculations
pdb_id="1a7i_A"
echo $pdb_id
input_file=${THIRD_PARTY_TOOLS}/check_setup/${pdb_id}.pdb
output_dir=${THIRD_PARTY_TOOLS}/check_setup/output/

. ${THIRD_PARTY_TOOLS}/runParty3Calc.sh $pdb_id $input_file $output_dir

