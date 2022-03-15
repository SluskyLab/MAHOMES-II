
import pandas as pd
import MAHOMES_II as MHMII
import sys

job_dir = str(sys.argv[1])
features = pd.read_csv("%s/features.csv"%(job_dir))
preds = MHMII.make_predictions_with_saved_MAHOMES_II(features)
preds.to_csv("%s/predictions.csv"%(job_dir))

