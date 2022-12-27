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
# @file   FeatureCalculations/features_pocket.py
# @brief  Combines all MAHOMES II pocket features (lining and void)
# @author Ryan Feehan (RFeehan93@gmail.com)

import pandas as pd
import numpy as np
import sys
import sys 

sys.path.insert(0, "%s" % "Pocket/")
import grid_tools as grid
import pocket_lining as lining
import features_ghecom as ghecom

def get_pocket_lining_features(this_data_dir, site, this_protein, res_nums):
    pdb = site.metal_atoms[0].struc_id
    metal_res=[]
    for metal in site.metal_atoms:
        metal_res.append(str(int(metal.seqID))+metal.chainID)

    SITE_center_pt = site.get_center()
    whole_pocket, pocket_info = ghecom.get_closest_pocket_coords(this_data_dir, pdb, SITE_center_pt)
    pocket_coords = whole_pocket[['X','Y','Z']]
    
    if len(whole_pocket) > 0:
        pocket_lining_res = lining.id_pocket_lining(pocket_coords, this_protein, cutoff=3)
        pocket_lining_res = set(pocket_lining_res).difference(metal_res)
        #print(pocket_lining_res, res_nums)
        
        pocket_lining_res = [x for x in pocket_lining_res if this_protein.residues[ this_protein.res_nums.index(x) ].type == "protein"]
        bb, bb_names, sc, sc_names = lining.split_pocket_lining(pocket_lining_res, pocket_coords, this_protein, cutoff = 2.2)
        #print(bb, bb_names, sc, sc_names)         
    else:
        print("No pocket", pdb, site.SITE_ID)
        bb_names = []
        sc_names = []
    
    labels, pocket_lining = lining.calc_lining_features(bb_names, sc_names, this_protein)
    pocket_lining = pd.DataFrame(pocket_lining, index=labels) 
    
    center_of_mass = np.mean(this_protein.Coords, axis = 0)
    dist_to_center = np.linalg.norm(SITE_center_pt - center_of_mass)
    labels, pocket_shape = grid.calc_pocket_features_ghecom(pocket_coords, SITE_center_pt, metal_res, this_protein, res_nums)
    labels.extend(["SITEDistCenter", "SITEDistNormCenter"])
    max_dist_to_center = np.sqrt(np.max( np.sum( (this_protein.Coords - dist_to_center)**2, axis = 1) ))
    pocket_shape.extend([dist_to_center, dist_to_center/max_dist_to_center])   
    pocket_shape = pd.DataFrame(pocket_shape, index=labels)

    ## add lining features that require Vol from pocket shape
    pocket_lining.loc["NoSC_vol"] = pocket_shape.loc["Vol"][0] + pocket_lining.loc["occ_vol"][0] 
    pocket_lining.loc["SC_vol_perc"] = pocket_lining.loc["occ_vol"][0] / pocket_lining.loc["NoSC_vol"][0]
    
    pocket_info = ghecom.calc_pocket_info_plus(pocket_coords, SITE_center_pt, center_of_mass, pocket_info)
    return(pocket_shape, pocket_lining, pocket_info)

## if ghecom can't find any pockets, the above willl fail and all (?) pocket features should be 0 in such an event
def get_empty_pocket_lining_features():
    pocket_shape_cols = ['Depth', 'Vol', 'LongPath', "SITEDistCenter", "SITEDistNormCenter",
    'farPtLow', 'PocketAreaLow', 'OffsetLow','LongAxLow', 'ShortAxLow',
    'farPtMid', 'PocketAreaMid','OffsetMid', 'LongAxMid', 'ShortAxMid',
    'farPtHigh','PocketAreaHigh', 'OffsetHigh', 'LongAxHigh', 'ShortAxHigh']
    pocket_lining_cols = ['in_pocket', 'num_pocket_bb', 'num_pocket_sc',
     'avg_eisen_hp','min_eisen', 'max_eisen', 'skew_eisen', 'std_dev_eisen',
     'avg_kyte_hp', 'min_kyte', 'max_kyte', 'skew_kyte','std_dev_kyte',
      'occ_vol', 'NoSC_vol', 'SC_vol_perc']
    pocket_info_cols = ['ClusterNum','Ngrid','Volume'
                   ,'Rinac(A) av','Rinac(A) mi','invRvolume(AA)'
                   ,'pocket_height','pocket_depth'
                   ,'metal_height','metal_depth'
                   ,'SITE_pocket_distance_min','num_adjacent_pockets'
                  ]
    pocket_shape = pd.DataFrame(0, columns=[0], index=pocket_shape_cols)
    pocket_lining = pd.DataFrame(0, columns=[0], index=pocket_lining_cols)
    pocket_info = pd.Series(0, index=pocket_info_cols)

    return(pocket_shape, pocket_lining, pocket_info)

def pocket_lining_features(struc_dir, site, this_protein, res_nums):
    #try:
    pocket_shape, pocket_lining, pocket_info = get_pocket_lining_features(struc_dir, site, this_protein, res_nums)
    #except:
    #    pocket_shape, pocket_lining, pocket_info = get_empty_pocket_lining_features()
    ## add prefixes for future domain handling and analysis
    pocket_shape = pocket_shape.T.add_prefix("pkt_shape_").T # double transform because rows can't prefix?
    pocket_lining = pocket_lining.T.add_prefix("pkt_lining_").T # double transform because rows can't prefix?
    pocket_info = pocket_info.add_prefix("pkt_info_")
    return(pocket_shape, pocket_lining, pocket_info)
