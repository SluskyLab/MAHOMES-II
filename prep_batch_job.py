# This script identifies all input structure files so they can e iterated over for feature calculations
#import pandas as pd
import glob
import os
import sys

job_dir = str(sys.argv[1])


with open("%s/batch_input.txt"%(job_dir), "w") as batch_input:
    for file_name in glob.glob("%s/*.pdb"%(job_dir)):
        pdb = file_name.split('/')[-1]
        pdb = pdb[:-4] # remove .pdb from file name
        #print(pdb)
        if len(pdb)>0:
            os.system("mkdir -p %s/%s"%(job_dir, pdb))
            batch_input.write("%s\n"%(pdb))

