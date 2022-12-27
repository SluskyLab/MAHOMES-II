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

import grid_tools
import numpy as np
import scipy
import warnings
from sklearn import cluster
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

eisenburg_hp = {"A":  0.620, "R": -2.530, "N": -0.780, "D": -0.900, "C":  0.290, "Q": -0.850, "E": -0.740, "G":  0.480, "H": -0.400,  "I":  1.380, 
"L":  1.060, "K": -1.500, "M":  0.640, "F":  1.190, "P":  0.120, "S": -0.180, "T": -0.050, "W":  0.810, "Y":  0.260, "V":  1.080, "X": 0}

kyte_doo_hp = {"A":  1.800, "R": -4.5, "N": -3.5, "D": -3.5, "C":  2.5, "Q": -3.5, "E": -3.5, "G":  -0.4, "H": -3.2,  "I":  4.5, 
"L":  3.8, "K": -3.9, "M":  1.9, "F":  2.8, "P":  -1.6, "S": -0.8, "T": -0.7, "W":  -0.9, "Y":  -1.2, "V":  4.2, "X": 0}

metal_size = {np.nan: 0,
         'MO':1, 'MOO':5, '4MO' : 1,'6MO' :1,'MOS': 4,
         'MG': 1,'3NI':1,'NI' : 1, 'ZN': 1,'MGF': 4,'MN3' : 1,'MN' : 1,'CO': 1,
         'OFE': 2, 'FE2': 1,'FEO': 3, 'FE' : 1,'FES': 4,
         'CU': 1, 'C2O' :3, 'CUA' : 2, 'CU1': 1,
         }

vdw_vols = {"G":48,"A":67,"S":73,"C":86,"P":90,"D":91,"T":93,"N":96,"V":105,
    "E":109,"Q":114,"H":118,"I":124,"L":124,"M":124,"K":135,"F":135,"Y":141,"R":148,"W":163,"X":0}
    

#added 12-18-18 MWF for clustering pocket points to identify the ones that are actually close to SITE center for identifying the residues lining the pocket
def cluster_grid_points(pocket_list):
    dbscan = cluster.DBSCAN(eps = 0.9)
    prediction = dbscan.fit_predict(pocket_list)
    return(prediction)

def graph_grid(pocket_list, prediction, near_neigh):
    colors = np.array(["black", "red", "green", "yellow", "orange", "purple"])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter3D(near_neigh[0], near_neigh[1], near_neigh[2], color="blue")
    ax.scatter3D(pocket_list[:,0], pocket_list[:,1], pocket_list[:,2], color=colors[prediction].tolist())
    plt.show()

def get_nearby_grid(pocket_list, cluster_prediction, SITE_center, cutoff):
    #cutoff should be in A
    grid_distances = np.sum((pocket_list - SITE_center)**2, axis=1)# calculate all distances between grid points and SITE center                
    metal_cluster = cluster_prediction[ grid_distances < float(cutoff)**2 ] #keep any cluster with a point within cutoff**2 of the site center (the distances are still squared)
    metal_cluster = sorted(set(metal_cluster))
    pocket_list = pocket_list[ np.in1d(cluster_prediction, metal_cluster) ]
    return(pocket_list)
    
def id_pocket_lining(grid_points, myProtein, cutoff=2.5):
    max_len=500
    if len(myProtein.residues) < max_len:
        all_dist = scipy.spatial.distance.cdist(myProtein.Coords, grid_points) #calculate all pairwise distances
        #print(np.amin(all_dist), np.amax(all_dist))
        pocket_liner = np.asarray(myProtein.Coords_index)[np.any(all_dist <= float(cutoff), axis = 1)].tolist() #keep the residues where at least one atom is within the cutoff distance
    # these PDBs are require to much memory, so they will be handled with a slower loop
    else:
        pocket_liner=[]
        for res in myProtein.residues:
            all_dist = scipy.spatial.distance.cdist(res.Coords, grid_points)        
            if np.any(all_dist <= float(cutoff)):
                res_num=str(res.resnum) + res.chain
                pocket_liner.append(res_num)

    pocket_liner = sorted(set(pocket_liner)) 
    #print(pocket_liner)
    return(pocket_liner)

def split_pocket_lining(pocket_lining_list, pocket_points, myProtein, cutoff = 2.2):
    backbone_res = []
    sidechain_res = []
    backbone_names = []
    sidechain_names = []
    #print(pocket_lining_list)
    #print(pocket_points)
    
    for entry in pocket_lining_list:
        this_res = myProtein.residues[ myProtein.res_nums.index(entry) ]
        #print(this_res.name, this_res.resnum)#, this_res.Atoms)
        if this_res.name in metal_size.keys():
            continue
        keep_atoms = this_res.Atoms.copy()
        if this_res.name != "GLY":
            keep_atoms = [x for x in keep_atoms if x not in ["H", "HA"]] #don't consider the hydrogens on the backbone
        this_dist = scipy.spatial.distance.cdist(this_res.Coords[ np.where(np.in1d(this_res.Atoms, keep_atoms)) ], pocket_points) #distance between all non-H coords of residue and the pocket grid points
        #print(this_dist.shape)
        if len(keep_atoms) < 5:
            continue
        closest_atoms = np.min(this_dist, axis = 1) #minimum distances for each atom
        min_bbone = np.min(closest_atoms[:4]) #minimum in first four atoms - C N CA O
        max_sc = np.min(closest_atoms[4:]) #minimum distance in sidechain
        if max_sc <= min_bbone: #sidechain is closer
            sidechain_res.append(entry) 
            sidechain_names.append(this_res.onelet)
            #print(this_res.onelet)
        elif max_sc < cutoff: #sidechain is pretty darn close anyways
            sidechain_res.append(entry)
            sidechain_names.append(this_res.onelet)
            #print(this_res.onelet)
        else: #otherwise, only backbone lines pocket
            backbone_res.append(entry)
            backbone_names.append(this_res.onelet)
            #print(this_res.onelet)
    #print(sidechain_names)
    #print(backbone_names)            
    return(backbone_res, backbone_names, sidechain_res, sidechain_names)

def calc_lining_features(bb_list, sc_list, myProtein):
    these_labels = ["in_pocket", "num_pocket_bb", "num_pocket_sc", "avg_eisen_hp", "min_eisen", "max_eisen", "skew_eisen", "std_dev_eisen", "avg_kyte_hp", "min_kyte", "max_kyte", "skew_kyte", "std_dev_kyte", "occ_vol"]
    if len(sc_list) > 0:
        eisenburg = np.asarray([eisenburg_hp[x] for x in sc_list])
        kyte = np.asarray([kyte_doo_hp[x] for x in sc_list])
        occ_vol = np.asarray([vdw_vols[x] for x in sc_list])
        eisen_descrip = scipy.stats.describe(eisenburg)
        kyte_descrip = scipy.stats.describe(kyte)
        these_features = [1, (len(bb_list) + len(sc_list)), len(sc_list), eisen_descrip.mean, eisen_descrip.minmax[0], eisen_descrip.minmax[1], eisen_descrip.skewness, np.sqrt(eisen_descrip.variance), 
            kyte_descrip.mean, kyte_descrip.minmax[0], kyte_descrip.minmax[1], kyte_descrip.skewness, np.sqrt(kyte_descrip.variance), sum(occ_vol)]
    else:
        these_features = [[0] * len(these_labels)]
    #print(these_features)
    return(these_labels, these_features)
    
    
