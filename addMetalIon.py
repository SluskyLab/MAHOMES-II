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

# @file  addMetalIon.py
# @brief The places a given metal ion average location of metal binding sidechain 
#        atom(s) (O, N, or S) for input coordinating residue ids and .pdb format structure.
#        This script was used to add metal ion binding sites for the AlphaFold set used in
#        MAHOMES II's publication (Feehan et al. 2023). Not all amino acid types are coded.
# @author Ryan Feehan <RFeehan93@gmail.com>

import glob
import os
import sys

import pandas as pd
import math
import gzip


## adjust pandas options so output pdb lines are missing anything
pd.set_option('display.max_columns',None)
pd.set_option('display.max_rows',None)
pd.set_option('display.width',None)
pd.set_option('display.max_colwidth',None)



####################################
#######   .pdb file parser   #######
####################################

## 'ATOM' and 'HETATM' use get_ATOM_DF for a atom level dataframe
def get_ATOM_DF(ATOM, pdb=None):
    atom_data = split_ATOM(ATOM['raw'])
    ATOM = pd.DataFrame({
         'serial':atom_data[0]  , 'name':atom_data[1]      ,  'altLoc':atom_data[2]     \
        ,'resName':atom_data[3] , 'chainID':atom_data[4]   ,  'seqNum':atom_data[5]     \
        ,'iCode':atom_data[6]   , 'x':atom_data[7]         ,  'y':atom_data[8]          \
        ,'z':atom_data[9]       , 'occupancy':atom_data[10] ,  'pLDDT':atom_data[11] \
        ,'segID':atom_data[12]   , 'element':atom_data[13]   ,  'charge':atom_data[14]     })
    
    if pdb!=None:
        ## add mapping variables
        ATOM['RES_ID'] = pdb[:4] + "_" + ATOM['chainID'].str[:] + "_" +  ATOM['resName'].str[:] + "_" + ATOM['seqNum'].str[:]
        ATOM['ATOM_ID'] =  ATOM['RES_ID'].str[:] + "_" + ATOM['element'].str[:]+ "_" + ATOM['serial'].str[:] 
    return(ATOM)


def split_ATOM(raw):
    return(  raw.str[6:11].str.replace(' ', ''),  raw.str[12:15].str.replace(' ', ''),     raw.str[16].str.replace(' ', ''), \
            raw.str[17:20].str.replace(' ', ''),     raw.str[21].str.replace(' ', ''),  raw.str[22:26].str.replace(' ', ''), \
               raw.str[27].str.replace(' ', ''),  raw.str[30:37].str.replace(' ', ''),  raw.str[38:45].str.replace(' ', ''), \
            raw.str[46:53].str.replace(' ', ''),  raw.str[54:59].str.replace(' ', ''),  raw.str[60:65].str.replace(' ', ''), \
            raw.str[72:75].str.replace(' ', ''),  raw.str[76:78].str.replace(' ', ''),  raw.str[79:].str.replace(' ', ''))

def get_AF_atom_df(AF_file):
    AF_file_df=pd.read_csv(AF_file, header=None, names = ['raw'])
    AF_file_df['raw'] = AF_file_df['raw'].str[:80]
    AF_file_df['key'] = AF_file_df['raw'].str[:6].str.replace(' ', '')
    AF_atom_df = get_ATOM_DF(AF_file_df.loc[AF_file_df['key']=="ATOM"])
    AF_atom_df.reset_index(drop=True, inplace=True)
    return(AF_atom_df)


#distance beween points, atoms in this case
def get_dis(x1, y1, z1, x2, y2, z2):
    X = float(x2) - float(x1)
    Y = float(y2) - float(y1)
    Z = float(z2) - float(z1)
    
    distance = math.sqrt((X*X) + (Y*Y) + (Z*Z))
    return(distance)

## returns the atomic coordinates most likely to be coordinating a metal ion 
## or the next best option depending on amino acid type
def get_resi_coords(cur_resi_df):
    new_atom_df = pd.DataFrame()
    if cur_resi_df['resName'].iloc[0]=="HIS":
        new_atom_df = cur_resi_df.loc[cur_resi_df['name'].isin(["NE"])]
    elif cur_resi_df['resName'].iloc[0]=="CYS":
        new_atom_df =cur_resi_df.loc[cur_resi_df['name'].isin(["SG"])]
    elif cur_resi_df['resName'].iloc[0]=="MET": 
        new_atom_df =cur_resi_df.loc[cur_resi_df['name'].isin(["SD"])]
    elif cur_resi_df['resName'].iloc[0] in ["THR", "SER"]:
        new_atom_df =cur_resi_df.loc[cur_resi_df['name'].isin(["OG"])]
    ## average coords of multiple atoms used for the following residues
    elif cur_resi_df['resName'].iloc[0] in ["GLU", "ASP", 'GLN', 'ASN', 'LEU', "PHE", "GLY"]:
        if cur_resi_df['resName'].iloc[0]=="GLU":
            cur_resi_df = cur_resi_df.loc[cur_resi_df['name'].isin(["OE"])]
        elif cur_resi_df['resName'].iloc[0]=="ASP":
            cur_resi_df = cur_resi_df.loc[cur_resi_df['name'].isin(["OD"])]
        elif cur_resi_df['resName'].iloc[0]=="ASN":
            cur_resi_df = cur_resi_df.loc[cur_resi_df['name'].isin(["ND","OD"])]
        elif cur_resi_df['resName'].iloc[0]=="GLN":
            cur_resi_df = cur_resi_df.loc[cur_resi_df['name'].isin(["NE","OE"])]
        ## just take the centroid of the ring for FEE and side chain C for LEU
        elif cur_resi_df['resName'].iloc[0] in ["PHE","LEU"]:
            cur_resi_df = cur_resi_df.loc[cur_resi_df['name'].isin(["CZ", "CE", "CD", "CG"])]
        ## use backbone N and O for GLY, which probably is more often correct for PHE and LEU too
        elif cur_resi_df['resName'].iloc[0]=="GLY":
            cur_resi_df = cur_resi_df.loc[cur_resi_df['name'].isin(["N","O"])]
        new_atom_df = cur_resi_df.iloc[:1].copy()
        new_atom_df['x'].iloc[:1] = cur_resi_df['x'].astype(float).mean()
        new_atom_df['y'].iloc[:1] = cur_resi_df['y'].astype(float).mean()
        new_atom_df['z'].iloc[:1] = cur_resi_df['z'].astype(float).mean()
    else:   
        print("amino acid type has not been coded yet")
        #print(cur_resi_df)
    return(new_atom_df)
            


def add_metal_to_AF(coord_resis, metal, input_file, output_file, chainID='A'):
    if len(coord_resis)==0:
        return(False)
    atom_df = get_AF_atom_df(input_file)    
    coord_atoms = []
    
    ## get coordinates to use for each coordinating resi
    for cur_resi in coord_resis:
        cur_resi_df = atom_df.loc[atom_df['seqNum']==cur_resi].copy()
        new_atom_df = get_resi_coords(cur_resi_df)
        if len(new_atom_df)>0:
            coord_atoms.append(new_atom_df)
    # finish making df with workable coordinates of coord resis 
    coord_resi_df = pd.concat(coord_atoms)
    for axis in ['x', 'y', 'z']: coord_resi_df[axis] = coord_resi_df[axis].astype(float)
    
    ## end without adding if not all expected resis were used or for the same metal ion site
    if len(coord_resi_df)!=len(coord_resis):
        print("a coordinating atom is missing for one of the given residues")
        return(False)
    for index1, row1 in coord_resi_df.iterrows():
        for index2, row2 in coord_resi_df.iterrows():
            if index1!=index2:
                distance = get_dis(row1['x'], row1['y'], row1['z'],row2['x'], row2['y'], row2['z'])
                if distance>12:
                    print("coordinating residues were not for same metal")
                    return(False)
    
    # calculate coordinates for metal
    metal_x = round(coord_resi_df['x'].astype(float).mean(),3)
    metal_y = round(coord_resi_df['y'].astype(float).mean(),3)
    metal_z = round(coord_resi_df['z'].astype(float).mean(),3)
    
    # find out what seqiID and serial should be used for the added metal ion
    metal_seqid = atom_df['seqNum'].astype(int).max()+1
    atom_num = atom_df['serial'].astype(int).max()+1
            
    # add metal to output file
    metal_line = "HETATM"+ str(atom_num).rjust(5)+" %s    %s "%(metal, metal)+ chainID + str(metal_seqid).rjust(4)+"    "+ str(metal_x).rjust(8)+ str(metal_y).rjust(8)+ str(metal_z).rjust(8)+"  1.00 99.99          %s\n"%(metal)
    metal_added=False
    with gzip.open(input_file, "rt") as orig_pdb:
        with open(output_file, "w+") as outData:
            for line in orig_pdb:
                outData.write(line)
                if line[0:3] == "TER":
                    if not metal_added:
                        outData.write(metal_line)
                        metal_added=True
                        
    return(metal_added)
