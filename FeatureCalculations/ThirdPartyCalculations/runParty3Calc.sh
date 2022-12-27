#!/bin/bash
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
#
# @file   FeatureCalculations/ThirdPartyCalculations/runParty3Calc.sh 
# @brief  Generates Rosetta, BLUUES, GHECOM, and FindGeo outputs required for MAHOMES II feature calculations
# @author Ryan Feehan <RFeehan93@gmail.com>
#
# usage: bash runParty3Calc.sh input_name path/to/input/file path/for/output path/to/rosetta path/to/bluues path/to/ghecom path/to/findgeo path/to/MAHOMES-II
# 
# Takes a .cif/.pdb file with a protein structure. Scores structure with Rosetta and outputs
# uniform .pdb file that is used for the rest of the feature calculations.
# Makes a .pqr file and uses it for bluues to make theoretical tirtration curve
# output. Outputs pocket grid(s) using ghecom. Then outputs coordination for 
# and metals in the file using findgeo.



echo "Working in runThirdPartyTools.sh!"

PDB=$1
INPUT_FILE=$2
OUTPUT_DIR=$3
echo $PDB
## get paths for required third paerty software
ROSETTA3=$4
BLUUES_DIR=$5
GHECOM_DIR=$6
MAHOMES_II_DIR=$7



## make sure we are using the correct environment
source ${MAHOMES_II_DIR}/venv/bin/activate

## prep and change to output directory
mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

THIRD_PARTY_TOOLS=${MAHOMES_II_DIR}/FeatureCalculations/ThirdPartyCalculations
#################################################################################
##################################   ROSETTA   ##################################
#################################################################################
## prep working directory
cp ${THIRD_PARTY_TOOLS}/Scoring.flags .
cp ${THIRD_PARTY_TOOLS}/EqualWeight.wts .

## for cluster
#module load rosetta/3.13
#rosetta_scripts.default.linuxgccrelease -s $INPUT_FILE -parser:protocol ${THIRD_PARTY_TOOLS}/ScoreOnly.xml @Scoring.flags > StdOutputScore.txt

## for server (and other machines)
${ROSETTA3}/bin/rosetta_scripts.default.linuxgccrelease -s $INPUT_FILE -parser:protocol ${THIRD_PARTY_TOOLS}/ScoreOnly.xml @Scoring.flags > StdOutputScore.txt

## delete the extra crap
mv ${PDB}_0001.pdb ${PDB}.pdb
rm Score.sc
rm Scoring.flags
rm EqualWeight.wts

#################################################################################
################################ bluues/pdb2pqr  ################################
#################################################################################
## prep a pdb file without any Rosetta scoring to avoid pdb2pqr errors
sed '/# All scores below are weighted scores/,$d' ${PDB}.pdb > ${PDB}_clean.pdb
python ${THIRD_PARTY_TOOLS}/PrepPDB_pdb2pqr.py ${OUTPUT_DIR}/${PDB}_clean.pdb
pdb2pqr30  --noopt --nodebump --ff=PARSE ${OUTPUT_DIR}/${PDB}_clean.pdb ${OUTPUT_DIR}/${PDB}.pqr --log-level CRITICAL
## remove cleaned pdb file
rm ${PDB}_clean.pdb

## run bluues 
${BLUUES_DIR}/bluues ${PDB}.pqr ${PDB}_bluues -pka > ${PDB}_bluues.log

#################################################################################
##################################   ghecom    ##################################
#################################################################################
${GHECOM_DIR}/ghecom -M M -ipdb ${PDB}.pdb -opocpdb ${PDB}_ghecom.pdb > ${PDB}_ghecom.log

#################################################################################
#################################   findgeo    ##################################
#################################################################################
## for cluster
#module load python/2.7
## otherwise
source ${MAHOMES_II_DIR}/venv2/bin/activate
## run findgeo with our hacks
python ${THIRD_PARTY_TOOLS}/findgeo.py -o -p ${PDB}.pdb -i ${PDB} -t 3.5 > ${PDB}_findgeo.log

