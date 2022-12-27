#!/bin/bash
#
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
# @file   FeatureCalculations/ThirdPartyCalculations/testRun.sh
# @brief  Makes third-party outputs for test pdb to compare with successful example outputs 
# @author Ryan Feehan <RFeehan93@gmail.com>

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

