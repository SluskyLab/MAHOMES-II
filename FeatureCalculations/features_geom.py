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
# @file   FeatureCalculations/features_geom.py
# @brief  Turns output from findgeo into MAHOMES II features
# @author Meghan W. Franklin
# @author Ryan Feehan (RFeehan93@gmail.com)

import os
import numpy as np
import pandas as pd


## this will all go in features_geom.py file
## after 38 in geom features have been changed on 8/30/2021 to remove charge, gRMSD, max gRMSD dev, valence, nVESCUM
# 0-35 geom classifiers, 36-38 irr/reg/distorted
# 39-42 N/O/S/other ligands,  43 overall RMSD
      
findgeo_geoms = ["lin", "trv", "tri", "tev", "spv", 
    "tet", "spl", "bva", "bvp", "pyv", 
    "spy", "tbp", "tpv", 
    "oct", "tpr", "pva", "pvp", "cof", "con", "ctf", "ctn",
    "pbp", "coc", "ctp", "hva", "hvp", "cuv", "sav",
    "hbp", "cub", "sqa", "boc", "bts", "btt", 
    "ttp", "csa", "irr"]


def collate_findgeo(these_metals, pdb_id, directory):
    #print(these_metals)
    geom_out_file = "%s/%s.findgeo"%(directory, pdb_id)
    this_geom = np.zeros(44)
    num_geoms = 0
    if os.path.isfile(geom_out_file) == True: 
        with open(geom_out_file, "r") as inData:
            #findgeo input is atomid, geom, irr/reg/distorted, RMSD, 4 N, O, S, other ligand ?
            for line in inData:
                line = line.strip().split("\t")
                ## findGeo can change FES to FE1 and FE2 and so forth, so do not include resName in check
                atom = line[0]
                if atom in these_metals:
                    num_geoms+=1
                    this_geom[ findgeo_geoms.index(line[1]) ] += 1 #this_geom[0:36] for the 36 regular geometries + irregular
                    if line[2] == "regular": 
                        this_geom[37] += 1
                        this_geom[43] += float(line[3]) #overal RMSD
                    elif line[2] == "distorted":
                        this_geom[38] += 1
                        this_geom[43] += float(line[3]) #overal RMSD
                    else:
                        # The findgeo.py script does not get the lowest RMSD when irr(!), but puts 100 instead so ... 1.5 might be close to max
                        this_geom[43] += 1.5 
                    this_geom[39] += float(line[4]) #N
                    this_geom[40] += float(line[5]) #O
                    this_geom[41] += float(line[6]) #S
                    this_geom[42] += float(line[7]) #other liganding atom
                    #print(float(line[3]), this_geom[43])
                    
    else:
        print("No geom data", pdb_id)

    these_labels = ["geom_" + x for x in findgeo_geoms]
    these_labels.extend(["geom_Reg", "geom_Distort", "geom_LigN", "geom_LigO", "geom_LigS", "geom_LigOther", "geom_AtomRMSD"])
    this_geom = pd.DataFrame.from_records(this_geom.reshape(-1, len(this_geom)).T, index = these_labels)
    
    this_geom = this_geom.T
    this_geom["geom_cn2"] = this_geom[["geom_lin", "geom_trv"]].sum(axis = 1)
    this_geom["geom_cn3"] = this_geom[["geom_tri", "geom_tev", "geom_spv"]].sum(axis = 1)
    this_geom["geom_cn4"] = this_geom[["geom_tet", "geom_spl", "geom_bva", "geom_bvp", "geom_pyv"]].sum(axis = 1)
    this_geom["geom_cn5"] = this_geom[["geom_spy", "geom_tbp", "geom_tpv"]].sum(axis = 1)
    this_geom["geom_cn6"] = this_geom[["geom_oct", "geom_tpr", "geom_pva", "geom_pvp", "geom_cof", "geom_con", "geom_ctf", "geom_ctn"]].sum(axis = 1)
    this_geom["geom_cn7"] = this_geom[["geom_pbp", "geom_coc", "geom_ctp", "geom_hva", "geom_hvp", "geom_cuv", "geom_sav"]].sum(axis = 1)
    this_geom["geom_cn8"] = this_geom[["geom_hbp", "geom_cub", "geom_sqa", "geom_boc", "geom_bts", "geom_btt"]].sum(axis = 1)
    this_geom["geom_cn9"] = this_geom[["geom_ttp", "geom_csa"]].sum(axis = 1)
    this_geom["geom_Filled"] = this_geom[["geom_lin", "geom_tri","geom_tet", "geom_spl", "geom_spy", "geom_tbp","geom_oct", "geom_tpr","geom_pbp", "geom_coc", "geom_ctp","geom_hbp", "geom_cub", "geom_sqa", "geom_boc", "geom_bts", "geom_btt","geom_ttp", "geom_csa"]].sum(axis = 1)
    this_geom["geom_PartFilled"] = this_geom[["geom_trv", "geom_tev", "geom_spv", "geom_bva", "geom_bvp", "geom_pyv", "geom_tpv", "geom_pva", "geom_pvp", "geom_cof", "geom_con", "geom_ctf", "geom_ctn", "geom_hva", "geom_hvp", "geom_cuv", "geom_sav"]].sum(axis = 1)
    
    if num_geoms > 0:
        this_geom["geom_AvgN"] = this_geom["geom_LigN"] / num_geoms
        this_geom["geom_AvgO"] = this_geom["geom_LigO"] / num_geoms
        this_geom["geom_AvgS"] = this_geom["geom_LigS"] / num_geoms
        this_geom["geom_AvgOther"] = this_geom["geom_LigOther"] / num_geoms
    else:
        this_geom["geom_AvgN"] = 0
        this_geom["geom_AvgO"] = 0
        this_geom["geom_AvgS"] = 0
        this_geom["geom_AvgOther"] = 0
    return(this_geom.T)
    
