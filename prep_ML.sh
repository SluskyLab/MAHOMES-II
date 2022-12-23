#!/bin/bash

## pass in local path to MAHOMES-II repo
MAHOMES_II_DIR=$1

echo "Running training and saving MAHOMES ML models with evaluations"
source ${MAHOMES_II_DIR}/venv/bin/activate
cd ${MAHOMES_II_DIR}/MachineLearning
python save_MAHOMES_II.py > saving_models_output.log

