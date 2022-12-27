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

## Rosetta does not output the serial for SSBOND, which is needed for pdb2pqr
import pandas as pd
import sys

pd.set_option('display.max_colwidth', 100)
pd.set_option('display.max_rows', 10000)

#####################################################################
######             CREATE WORKABLE PDB DATAFRAME               ######
#####################################################################
PDB_file = str(sys.argv[1])

## simple dataframe from pdb (with NMR adjustment)
def fix_SSBOND(PDB_file):
    f = open(PDB_file, "r")
    pdb_lines = f.readlines()
    pdb_file_df = pd.DataFrame(pdb_lines, columns = ['raw'])
    pdb_file_df['raw'] = pdb_file_df['raw'].str[:80]
    pdb_file_df['key'] = pdb_file_df['raw'].str[:6].str.replace(' ', '')
    #print(len(pdb_file_df.loc[pdb_file_df['key'].str[:]=="SSBOND"]))
    if len(pdb_file_df.loc[pdb_file_df['key'].str[:]=="SSBOND"])<1:
        return
    elif len(pdb_file_df.loc[pdb_file_df['key'].str[:]=="SSBOND"])<100:
        i=0
        for index, row in pdb_file_df.loc[pdb_file_df['key'].str[:]=="SSBOND"].iterrows():
            i+=1
            new_line =row['raw'][:9]+str(i)+row['raw'][10:]
            if len(str(i))>1:
                new_line = row['raw'][:8]+str(i)+row['raw'][10:]
            #print(pdb_file_df.loc[index]['raw'])
            #new_line = row['raw'][:9]+str(i)+row['raw'][10:]
            #print(new_line)
            pdb_file_df.loc[index]['raw']=new_line

        del  pdb_file_df['key']
        pdb_file_df.to_string(PDB_file, header=False, index=False, col_space=80)
fix_SSBOND(PDB_file)
