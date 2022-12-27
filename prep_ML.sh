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

# @file  
# @brief 
# @author Ryan Feehan <RFeehan93@gmail.com>

## pass in local path to MAHOMES-II repo
MAHOMES_II_DIR=$1

echo "Running training and saving MAHOMES ML models with evaluations"
source ${MAHOMES_II_DIR}/venv/bin/activate
cd ${MAHOMES_II_DIR}/MachineLearning
python save_MAHOMES_II.py > saving_models_output.log

