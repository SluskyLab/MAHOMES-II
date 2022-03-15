## Rosetta does not output the serial for SSBOND, which is needed for pdb2pqr
import pandas as pd
import sys

pd.set_option('max_colwidth', 100)
pd.set_option('max_rows', 10000)

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
    print(len(pdb_file_df.loc[pdb_file_df['key'].str[:]=="SSBOND"]))
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
