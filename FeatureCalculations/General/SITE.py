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

import numpy as np
import pandas as pd
import math
import glob
import os
import sys

import PDBparser_dataframe as PDB

#pdb = str(sys.argv[1])

## for identifying ions of interest
PDB_ION_CODES = [ "FE" ,"FE2","FES","FEO"
                 ,"CU" ,"CU1","CUA"
                 ,"MG" ,"ZN" ,"MN"
                 ,"MO" ,"MOO","MOS"
                 ,"NI" ,"3CO","CO"]

## for identifying atoms of interest
ION_PERIODIC_LIST = ["FE", "CU", "MG", "ZN", "MN", "MO", "NI", "CO"]


#distance beween points, atoms in this case
def get_dis(x1, y1, z1, x2, y2, z2):
    X = float(x2) - float(x1)
    Y = float(y2) - float(y1)
    Z = float(z2) - float(z1)
    
    distance = math.sqrt((X*X) + (Y*Y) + (Z*Z))
    return(distance)

class ATOM:
    def __init__(self, struc_id, chainID, seq_num, res_name, atom_serial, atom_name, x, y, z, site_id=None):
        self.struc_id = struc_id
        self.chainID = chainID
        self.seqID = seq_num
        self.resName = res_name
        self.serial = atom_serial
        self.atom_name = atom_name
        self.X = float(x)
        self.Y = float(y)
        self.Z = float(z)
        site_id = site_id
        
    
    def distance_to_COORD(self, x, y, z):
        distance = get_dis(self.X, self.Y, self.Z, x, y, z)
        return(distance)
    
    def distance_to_ATOM(self, atom):
        distance = self.distance_to_COORD(atom.X, atom.Y, atom.Z)
        return(distance)

class SITE:
    def __init__(self, site_id=None):
        self.metal_atoms = []
        self.SITE_ID = site_id
        self.bad=False
        self.note=""
        
    def get_struc_id(self):
        if len(self.metal_atoms)>0:
            return(self.metal_atoms[0].struc_id)
        else:
            return("")

    def add_metal_ATOM(self, metal_atom):
        metal_atom.SITE_ID = self.SITE_ID
        self.metal_atoms.append(metal_atom)
    
    
    def distance_to_COORD(self, x, y, z):
        distance = self.metal_atoms[0].distance_to_COORD(x, y, z)
        if len(self.metal_atoms) == 1:
            return(distance)
        
        for metal_atom in self.metal_atoms[1:]:
            new_dist = metal_atom.distance_to_COORD(x, y, z)
            if distance >= new_dist:
                distance=new_dist
        return(distance)
    
    def distance_to_ATOM(self, atom):
        distance = self.distance_to_COORD(atom.X, atom.Y, atom.Z)
        return(distance)
    
    def distance_to_SITE(self, site):
        distance = self.distance_to_ATOM(site.metal_atoms[0])
        if len(site.metal_atoms) == 1:
            return(distance)
        
        for metal_atom in site.metal_atoms[1:]:
            new_dist = self.distance_to_ATOM(metal_atom)
            if distance >= new_dist:
                distance=new_dist
        return(distance)

    def get_center(self):
        SITE_X=0; SITE_Y=0; SITE_Z =0
        for metal in self.metal_atoms:
            SITE_X+=metal.X; SITE_Y+=metal.Y; SITE_Z+=metal.Z
        SITE_X=SITE_X/len(self.metal_atoms)
        SITE_Y=SITE_Y/len(self.metal_atoms)
        SITE_Z=SITE_Z/len(self.metal_atoms)
        return([SITE_X, SITE_Y, SITE_Z])

    def check_site_quality(self, atom_df):
        close_cutoff=6
        ## check for any close residues
        SITE_center_pt = self.get_center()
        atom_df['dist'] = np.sqrt( (atom_df.x-SITE_center_pt[0])**2 + (atom_df.y-SITE_center_pt[1])**2 + (atom_df.z-SITE_center_pt[2])**2)
        close_atoms_df = atom_df.loc[atom_df['dist']<=close_cutoff]
        if 1>len(close_atoms_df):
            self.bad=True
            self.note="no residues within %s angstroms of site center"%str(close_cutoff)

    def get_residue_distances(self, protein):
        distances = []
        for res in range(0, len(protein.residues)):
            this_res = protein.residues[res]
            if "protein"==this_res.type:
                res_coords = [0, 0, 0]; num_atoms=0
                resX=0;resY=0;resZ=0
                for atom in this_res.Atoms:
                    if atom in ['CA', 'CB', 'CD', 'CD1','CD2','CE','CE1','CE2','CE3','CG','CG1','CG2','CH2','CZ','CZ2','CZ3']:
                        num_atoms+=1
                        resX+=this_res.Coords[ this_res.Atoms.index(atom) ][0]
                        resY+=this_res.Coords[ this_res.Atoms.index(atom) ][1]
                        resZ+=this_res.Coords[ this_res.Atoms.index(atom) ][2]
                
                if num_atoms ==0:
                    res_coords=[0,0,0]
                    if len(this_res.Coords)==1:
                        res_coords=this_res.Coords[0]
                    else:
                        res_coords=np.mean(this_res.Coords)
                    resX = res_coords[0]; resY = res_coords[1]; resZ = res_coords[2]
                else:
                    resX=resX/num_atoms; resY=resY/num_atoms; resZ=resZ/num_atoms
                this_distance = self.distance_to_COORD(float(resX), float(resY), float(resZ))
                distances.append(pd.DataFrame(
                    [[this_res.name+'_'+str(this_res.res_index), str(this_res.resnum), this_res.chain, this_distance]]
                    ,columns=["ResID", "seqID", "chainID", "C_Dist"]
                ))
        distances = pd.concat(distances, ignore_index=True)
        return(distances)


##########################################
## find all sites for a given structure ##
##########################################

def reformat_output_to_SITE(dataset, struc_id):
    input_dataset = dataset.copy()
    all_sites = []
    for site_id in input_dataset['SITE_ID'].unique().tolist():
        new_site = SITE(site_id)
        site_set = input_dataset[input_dataset['SITE_ID']==site_id].copy()
        site_set=site_set.sort_values(['seqNum', 'serial'], ascending=[True, True]).reset_index(drop=True)
        for index, row in site_set.iterrows():
            new_atom = ATOM(struc_id, row.chainID, row.seqNum, row.resName, row.serial, row.atom_name, row.x, row.y, row.z)
            new_site.add_metal_ATOM(new_atom)
        all_sites.append(new_site)
    return(all_sites)

def find_SITEs(pdb_id, input_file):
    ## create pdb ligand dataframe 
    pdb_df= PDB.get_PDB_DF(pdb_id, input_file)
    pdb_df = PDB.get_MODRES_pdb_df(pdb_df)
    hetatm_df = PDB.get_ATOM_DF(pdb_df[pdb_df['key']=='HETATM'], pdb_id)
    ## get rid of everything but metal atoms/ions of interest
    hetatm_df = hetatm_df[hetatm_df['element'].isin(ION_PERIODIC_LIST)]
    hetatm_df = hetatm_df[hetatm_df['resName'].isin(PDB_ION_CODES)]
    hetatm_df = hetatm_df.reset_index(drop=True)
    
    ############################################
    ## assign SITE numbers to metal ion atoms ##
    ############################################
    hetatm_df['site_num'] = -1
    for curAtom_i, curAtom_row in hetatm_df.iterrows():
        curAtom = ATOM(pdb_id, curAtom_row.chainID, curAtom_row.seqNum, curAtom_row.resName, curAtom_row.serial, curAtom_row.atom_name, curAtom_row.x, curAtom_row.y, curAtom_row.z)
        ## assign atom to new SITE num if it does not already have a SITE num
        curAtom_site_num = hetatm_df.loc[curAtom_i].site_num
        if curAtom_site_num == -1:
            curAtom_site_num = hetatm_df['site_num'].max()+1
            hetatm_df.loc[curAtom_i, 'site_num'] = curAtom_site_num
        
        ## check for other atoms that are close to be added to same SITE
        for newAtom_i, newAtom_row in hetatm_df[curAtom_i+1:].iterrows():
            newAtom_distance = curAtom.distance_to_COORD(newAtom_row.x, newAtom_row.y, newAtom_row.z)
            if newAtom_distance <= 5.0:
                newAtom_site_num = hetatm_df.loc[newAtom_i].site_num
                if (newAtom_site_num == -1):
                    ## newAtom has not yet been assigned to SITE
                    hetatm_df.loc[newAtom_i, 'site_num'] = curAtom_site_num
                else:
                    ## change all atoms added to curAtom's site_num to previously made site_num
                    hetatm_df.loc[hetatm_df['site_num']==curAtom_site_num, 'site_num'] = newAtom_site_num
                    ## update cur_site_num so other ions will be added to proper site_num
                    curAtom_site_num = newAtom_site_num
                    
    hetatm_df['SITE_ID'] = pdb_id + "_" + hetatm_df['site_num'].astype(str)

    return_sites = reformat_output_to_SITE(hetatm_df, pdb_id)
    atom_df = PDB.get_ATOM_DF(pdb_df[pdb_df['key']=='ATOM'], pdb=None)
    for site in return_sites:
        site.check_site_quality(atom_df)
    return(return_sites)

def make_known_SITE(pdb_id, input_file, prev_site_id, site_metal_SeqIDs):
    ## create pdb ligand dataframe 
    pdb_df= PDB.get_PDB_DF(pdb_id, input_file)
    pdb_df = PDB.get_MODRES_pdb_df(pdb_df)
    hetatm_df = PDB.get_ATOM_DF(pdb_df[pdb_df['key']=='HETATM'], pdb_id)
    ## get rid of everything but metal atoms/ions of interest
    hetatm_df = hetatm_df[hetatm_df['element'].isin(ION_PERIODIC_LIST)]
    hetatm_df = hetatm_df[hetatm_df['resName'].isin(PDB_ION_CODES)]
    hetatm_df = hetatm_df.loc[hetatm_df['seqNum'].astype(int).isin(site_metal_SeqIDs)]
    hetatm_df = hetatm_df.reset_index(drop=True)
    
    ############################################
    ## assign SITE numbers to metal ion atoms ##
    ############################################
    hetatm_df['site_num'] = prev_site_id[5:]
    hetatm_df['SITE_ID'] = pdb_id + "_" + hetatm_df['site_num'].astype(str)

    return_sites = reformat_output_to_SITE(hetatm_df, pdb_id)
    atom_df = PDB.get_ATOM_DF(pdb_df[pdb_df['key']=='ATOM'], pdb=None)
    for site in return_sites:
        site.check_site_quality(atom_df)
    return(return_sites)

