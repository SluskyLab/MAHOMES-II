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

# @file   prep_batch_job.py
# @brief  Identifies structure files located in the given input directory
# @author Ryan Feehan <RFeehan93@gmail.com>

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

    for file_name in glob.glob("%s/*.cif"%(job_dir)):
        pdb = file_name.split('/')[-1]
        pdb = pdb[:-4] # remove .mmcif from file name
        if len(pdb)>0:
            os.system("mkdir -p %s/%s"%(job_dir, pdb))
            batch_input.write("%s\n"%(pdb))

