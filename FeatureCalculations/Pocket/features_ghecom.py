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
# @file   FeatureCalculations/Pocket/features_ghecom.py 
# @brief  Finds site's pocket grid from GHECOM output and adds some MAHOMES II pocket void features
# @author Ryan Feehan (RFeehan93@gmail.com)

import pandas as pd
import numpy as np
import math
import sys
import os
import sys 
import scipy
import subprocess

RESOURCE_DIR = "P/"
sys.path.insert(0, "%s" % "Pocket/")
import grid_tools as grid
import pocket_lining as lining

def get_ghecom_pocket_info(filename):
    cluster_cens = subprocess.check_output(["grep", "CA  CEN", filename])
    cluster_cens = cluster_cens.decode("utf-8").strip().split("\n")
    new_cens = []
    for line in cluster_cens:
        new_cens.append( [ line[20:22].strip(), float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())] )
    cluster_cens = pd.DataFrame.from_records(new_cens, columns = ["PocketID", "X", "Y", "Z"])
    num_clusters = len(cluster_cens)
    cluster_cens["ClusterNum"] = ""
    cluster_cens["Ngrid"] = 0
    cluster_cens["Volume"] = 0
    cluster_cens["Rinac(A) av"] = 0
    cluster_cens["Rinac(A) mi"] = 0
    cluster_cens["invRvolume(AA)"] = 0
    if len(cluster_cens)==0:
        return(cluster_cens)
    #charmingly, the chain IDs that ghecom outputs don't match the cluster numbers after 9 - the chains start at A and the clusters numbers continue. Wheeeee.
    keep_clusters = cluster_cens.PocketID.values
    good_clusters = []
    for x in keep_clusters:
        try:
            good_clusters.append(str(int(x)))
            cluster_cens.loc[cluster_cens["PocketID"].astype(str) == str(x), "ClusterNum"] = str(int(x))
        except ValueError:
            good_clusters.append(str(ord(x)-55)) #convert the letters to numbers
            cluster_cens.loc[cluster_cens["PocketID"].astype(str) == str(x), "ClusterNum"] = str(ord(x)-55)

    small_clusters = subprocess.check_output( ["grep", "REMARK  CLUSTER_PROPERTY  ", filename ])
    small_clusters = small_clusters.decode("utf-8").strip().split("\n")
    small_clusters = [ x.split() for x in small_clusters[0:num_clusters] ]
    for x in small_clusters:
        #print(x)
        if x[2] in cluster_cens.ClusterNum.values:
            cluster_cens.loc[cluster_cens["ClusterNum"] == x[2], "Ngrid"] = float(x[4])
            cluster_cens.loc[cluster_cens["ClusterNum"] == x[2], "Volume"] = float(x[6])
            cluster_cens.loc[cluster_cens["ClusterNum"] == x[2], "Rinac(A) av"] = float(x[9])
            cluster_cens.loc[cluster_cens["ClusterNum"] == x[2], "Rinac(A) mi"] = float(x[11])
            cluster_cens.loc[cluster_cens["ClusterNum"] == x[2], "invRvolume(AA)"] = float(x[13])
    return(cluster_cens)


def get_clusters_COORDs(grid_file):
    read_data=True
    het_lines = []
    with open(grid_file, "r") as inData:
        for line in inData:
            if line[0:3] == "TER":
                read_data = False
            elif line[0:17] == "REMARK  CLUSTER  ":
                read_data=True
            if read_data == True:
                if line[0:6] == "HETATM":
                    if line[12:16].strip() != "CA":
                        line = line.strip()
                        het_lines.append([ line[20:22].strip(), float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ])
    grid_coords = pd.DataFrame.from_records(het_lines, columns = ["PocketID", "X", "Y", "Z"])
    return(grid_coords)
                        
def get_closest_pocket_coords(directory, pdb, SITE_center_pt):
    ## read in information from ghecom output for this structure
    ghecom_output_f = "%s/%s_ghecom.pdb"%(directory, pdb)
    
    grid_coords = get_clusters_COORDs(ghecom_output_f)
    pockets_info = get_ghecom_pocket_info(ghecom_output_f)
    
    ## find closest pocket to metal(s)
    grid_coords['dist'] = np.sqrt( (grid_coords.X-SITE_center_pt[0])**2 + (grid_coords.Y-SITE_center_pt[1])**2 + (grid_coords.Z-SITE_center_pt[2])**2)
    closest_dist= grid_coords['dist'].min()
    closest_pocket = grid_coords.loc[grid_coords['dist']==closest_dist,'PocketID'].values[0]
    
    ## find the (maybe) better pocket when multiple are adjacent to site
    adjacent_pockets = grid_coords.loc[grid_coords['dist']<=5,'PocketID'].unique()
    if 1<len(adjacent_pockets):
        ## get the closest "large" pocket when possable
        tmp_pocket_info = pockets_info.loc[pockets_info['PocketID'].isin(adjacent_pockets)].copy()
        tmp_pocket_info = tmp_pocket_info.loc[tmp_pocket_info['Volume']>100]
        if len(tmp_pocket_info['PocketID'].unique())>0:
            grid_coords = grid_coords.loc[grid_coords['PocketID'].isin(tmp_pocket_info['PocketID'].unique())]
            closest_dist= grid_coords['dist'].min()
            closest_pocket = grid_coords.loc[grid_coords['dist']==closest_dist,'PocketID'].values[0]

        
    ## return grid and ghecom info for closest pocket
    return_grid_coords = grid_coords.loc[grid_coords['PocketID']==closest_pocket].copy()
    pocket_info = pockets_info.loc[pockets_info['PocketID']==closest_pocket].copy()
    pocket_info['SITE_pocket_distance_min']=closest_dist
    pocket_info['num_adjacent_pockets']=len(adjacent_pockets)
    return(return_grid_coords.reset_index(drop=True), pocket_info.reset_index(drop=True).iloc[0])


def rotate_pocket_grid(coords_df, centroid1, centroid2):
    #rotate and translate coords to the new centroid1-centroid2 axis
    axis_length = np.linalg.norm(np.asarray(centroid1) - np.asarray(centroid2))
    new_coords = grid.rotate_trans_coords(np.asarray(centroid1), np.asarray(centroid2), np.asarray([0,0,0]), np.asarray([0,0,axis_length]), coords_df[["X", "Y", "Z"]])
    new_coords = pd.DataFrame(np.asarray(new_coords), columns=['X', 'Y', 'Z'])
    return(new_coords)

def calc_distance(p1, p2):
    return( math.sqrt( ((p1[0]-p2[0])*(p1[0]-p2[0])) + ((p1[1]-p2[1])*(p1[1]-p2[1])) + ((p1[2]-p2[2])*(p1[2]-p2[2])) ) )

def calculate_pocket_height_and_depth(coords_df, metal_coords):
    #calculate max distances; after rotation to a z-axis, this is the max - min of the z-direction
    max_z = coords_df.loc[coords_df['Z'].argmax()]
    min_z = coords_df.loc[coords_df['Z'].argmin()]
    #print(max_z)
    #print(min_z)
    height = max_z.Z-min_z.Z
    #print(height)
    depth = calc_distance(max_z, min_z)
    #print(depth)
    metal_hieght = calc_distance(metal_coords, min_z)
    metal_depth = calc_distance(metal_coords, max_z)
    return( height, depth, metal_hieght, metal_depth )


#added 11/3/2021 RF
def calc_pocket_info_plus(pocket, site_center, center_of_mass, pocket_info):
    ## rotate pocket so protein center is oriding and z axis points to pocket center
    pocket_center = np.mean(pocket[["X", "Y", "Z"]], axis = 0)
    pocket_w_site=pd.DataFrame([site_center], columns=['X', 'Y', 'Z'])
    pocket_w_site=pd.concat([pocket_w_site, pocket], ignore_index=True)
    new_pocket_w_site = rotate_pocket_grid(pocket_w_site, center_of_mass, pocket_center)
    new_site_center = new_pocket_w_site.iloc[0].copy()
    new_pocket = new_pocket_w_site.iloc[1:].copy()
    new_pocket_center = np.mean(new_pocket[["X", "Y", "Z"]], axis = 0) # get new pocket center for calculations
    pocket_info['pocket_height'], pocket_info['pocket_depth'], pocket_info['metal_height'], pocket_info['metal_depth'] = calculate_pocket_height_and_depth(new_pocket, new_site_center)
    del pocket_info['X'];del pocket_info['Y'];del pocket_info['Z']; del pocket_info['PocketID']
    return(pocket_info)

