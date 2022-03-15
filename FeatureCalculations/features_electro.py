## this will all go in features_electro.py file
import os
import sys
import numpy as np
import pandas as pd
import scipy



electro_inCutoff = 3.5
electro_outCutoff = 9


moment_terms=['mu2','mu3','mu4']
other_electro_terms = ["GBR6", "pKa_shift", "dpKa_titr",  "dpKa_desolv", "dpKa_bg", "DestabRank", "StabRank", "ResSigDev", "SolvEnergy", "SolvExp"]
electro_variables = moment_terms.copy()
for term in other_electro_terms: electro_variables.append(term)

def extract_bluues_data(df_in, prefix):
    df = df_in[electro_variables].copy()
    
    this_mean = df.mean(axis = 0, numeric_only = True)
    this_mean = this_mean.add_prefix("Elec_%s_mean_"%prefix)
    
    this_std = df.std(axis = 0, numeric_only = True)
    this_std.fillna(0, inplace=True) # will be nan for cases with 1 residue
    this_std = this_std.add_prefix("Elec_%s_std_"%prefix)
    
    this_max = df.max(axis = 0, numeric_only=True)
    this_max = this_max.add_prefix("Elec_%s_max_"%prefix)
    
    this_min = df.min(axis = 0, numeric_only = True)
    this_min = this_min.add_prefix("Elec_%s_min_"%prefix)

    this_range = df.max(axis = 0, numeric_only = True) - df.min(axis = 0, numeric_only = True)
    this_range = this_range.add_prefix("Elec_%s_range_"%prefix)
    
    bluues_data = pd.concat([this_mean, this_std, this_min, this_max, this_range])
    
    if 0==len(df_in):bluues_data.fillna(0, inplace=True) # no ion residues will make almost everything NaN
    for term in moment_terms:
        bluues_data['Elec_%s_Z_%s_count'%(prefix, term)]=len(df_in.loc[df_in['Z_%s'%term]>1])
        
    return(bluues_data)


def z_n(X, n):
    z_n = (X[n]-np.mean(X))/np.std(X)
    return(z_n)

def add_z_scores(titr_data):
    for term in moment_terms:
        titr_data['Z_%s'%term]=0
        this_mean = titr_data[term].mean()
        this_std = titr_data[term].std()
        for index, row in titr_data.iterrows():
            titr_data.loc[index, 'Z_%s'%term]=np.abs((row[term]-this_mean)/this_std)
    return(titr_data)

def get_xenvr(df_in):
    env_df = df_in.copy()
    env_df.loc[:,"dist_weight"] = 1/(env_df.loc[:,"C_Dist"]**2)
    
    for term in electro_variables:
        env_df[term] = env_df[term]*env_df["dist_weight"]
        
    env_df = env_df[electro_variables]
    xenvr = env_df.sum(axis = 0, numeric_only = True)
    xenvr = xenvr.add_prefix("Elec_xenv_")
    return(xenvr)
def get_titr_data(pdb_id, directory, distances):
    elec_out_file = "%s/%s_ElectFeatures.txt"%(directory, pdb_id)
    if os.path.isfile(elec_out_file) == True: #is there bluues electrostatics data? #possibly change this to calculate features if possible and np.nan if no data
        titr_data = pd.read_csv(elec_out_file, sep = "\t", header = 0)
    else:
        print("No bluues data", pdb_id)
        return(pd.DataFrame())
    num_titr_resis=len(titr_data)
    titr_data.rename(index = str, columns = {"Unnamed: 0":"ResID"}, inplace = True)
    titr_data.loc[:,'ResName'] = titr_data.loc[:, 'ResID'].str[:3]
    titr_data = add_z_scores(titr_data)
    titr_data.loc[:,"seqID"] = titr_data.loc[:,"ResID"].str[4:]
    distances = distances[["seqID", "C_Dist"]]
    titr_data = pd.merge(titr_data, distances, on = "seqID", how = "inner")
    if len(titr_data)<num_titr_resis:
        print("check %s for ionizable residues lost during electro calc. prep"%(pdb_id))
    return(titr_data)

def collate_bluues(pdb_id, directory, neighbor_list1, neighbor_list2, distances):
    titr_data = get_titr_data(pdb_id, directory, distances)
    neighbors = [x[:-1] for x in neighbor_list1] #drop chain from numbers
    inside_neigh = titr_data[titr_data.seqID.isin(neighbors)]

    neighbors2 = set(neighbor_list2).difference(neighbor_list1)
    neighbors2 = [x[:-1] for x in neighbors2]
    outside_neigh = titr_data[titr_data.seqID.isin(neighbors2)]

    inside = extract_bluues_data(inside_neigh, "ins")
    outside = extract_bluues_data(outside_neigh, "outs")
    
    neighbors_9a = [x[:-1] for x in neighbor_list2] #drop chain from numbers
    neighbors_9a_data = titr_data[titr_data.seqID.isin(neighbors_9a)]
    xenvr = get_xenvr(neighbors_9a_data)
    this_data = pd.concat([inside, outside, xenvr])
    this_data['Elec_ins_numResis']=len(inside_neigh)
    this_data['Elec_outs_numResis']=len(outside_neigh)
    this_data['Elec_xenvr_numResis']=len(neighbors_9a_data)
    #print(this_data)
    return(this_data)



def get_electro_features(directory, site, this_protein):
    pdb_id=site.metal_atoms[0].struc_id
    SITE_center_pt=site.get_center()
    metal_seqIDs=[]
    for metal in site.metal_atoms:
        metal_seqIDs.append(metal.seqID)

    #"%s/%s_ElectFeatures.txt"%(directory, pdb_id)
    #electrostatics calculations
    inner_shell = sorted(set(this_protein.get_neighbor_res(SITE_center_pt, electro_inCutoff)).difference(metal_seqIDs))
    outer_shell = sorted(set(this_protein.get_neighbor_res(SITE_center_pt, electro_outCutoff)).difference(metal_seqIDs))

    if len(metal_seqIDs)>1:
        inner_shell =  []; outer_shell=[]
        for metal in site.metal_atoms:
            metal_coords = [metal.X, metal.Y, metal.Z]
            new_inner = sorted(set(this_protein.get_neighbor_res(metal_coords, electro_inCutoff)).difference(metal_seqIDs))
            for r in new_inner: 
                if r not in inner_shell: inner_shell.append(r)
            new_outer = sorted(set(this_protein.get_neighbor_res(metal_coords, electro_outCutoff)).difference(metal_seqIDs))
            for r in new_outer: 
                if r not in outer_shell: outer_shell.append(r)
    distances = site.get_residue_distances(this_protein)
    bluues_electro = collate_bluues(pdb_id, directory, inner_shell, outer_shell, distances)
    return(bluues_electro)







