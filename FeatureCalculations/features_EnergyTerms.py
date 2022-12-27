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
# @file   FeatureCalculations/features_EnergyTerms.py
# @brief  Turns reformated output per residue Rosetta scores into MAHOMES II features
# @author Ryan Feehan (RFeehan93@gmail.com)


import numpy as np
import pandas as pd
import scipy

#scoring terms for labeling    
energy_cols= ['fa_atr', 'fa_rep', 'fa_sol', 'fa_intra_atr_xover4',
               'fa_intra_rep_xover4', 'fa_intra_sol_xover4', 'lk_ball', 'lk_ball_iso',
               'lk_ball_bridge', 'lk_ball_bridge_uncpl', 'fa_elec', 'fa_intra_elec',
               'pro_close', 'hbond_sr_bb', 'hbond_lr_bb', 'hbond_bb_sc', 'hbond_sc',
               'dslf_fa13', 'omega', 'fa_dun', 'fa_dun_dev', 'fa_dun_rot',
               'fa_dun_semi', 'p_aa_pp', 'yhh_planarity', 'hxl_tors', 'ref','rama_prepro'
               ,'total']

## returns dataframe with Rosetta residue scores 
def get_res_scores(directory, pdb_id):
    ## read in Rosetta score file
    filename1 = "%s/%s_ByResValues.txt"%(directory, pdb_id)
    by_res_scores =  pd.read_csv(filename1, sep='\t')

    seq_ids = by_res_scores["label"].values
    seq_ids = [x.split("_")[-1] for x in seq_ids]
    by_res_scores["ResNum"] = seq_ids
    by_res_scores["resName"]=by_res_scores["label"].str[:3]
    by_res_scores["ResID"]=by_res_scores["resName"]+"_"+by_res_scores["ResNum"]
    return(by_res_scores)
    
def calculate_dist_to_Cs(protein, atom_coords):
    atom_type_coords = np.zeros((len(protein.residues), 3))
    label_list = []
    
    for res in range(0, len(protein.residues)):
        this_res = protein.residues[res]
        new_label = this_res.name+"_"+str(this_res.res_index+1)+","+str(this_res.resnum)+","+str(this_res.type)
        label_list.append(new_label)
        
        res_coords = [0, 0, 0]; num_atoms=0
        for atom in this_res.Atoms:
            #print(type(atom))
            if atom in ['CA', 'CB', 'CD', 'CD1','CD2','CE','CE1','CE2','CE3','CG','CG1','CG2','CH2','CZ','CZ2','CZ3']:
                num_atoms+=1
                res_coords[0]+=this_res.Coords[ this_res.Atoms.index(atom) ][0]
                res_coords[1]+=this_res.Coords[ this_res.Atoms.index(atom) ][1]
                res_coords[2]+=this_res.Coords[ this_res.Atoms.index(atom) ][2]
        
        if num_atoms ==0:
            res_coords=np.mean(this_res.Coords)
        else:
            res_coords[0]=res_coords[0]/num_atoms
            res_coords[1]=res_coords[1]/num_atoms
            res_coords[2]=res_coords[2]/num_atoms
        atom_type_coords[res] = res_coords
        
    all_dist = scipy.spatial.distance.cdist(atom_type_coords, [atom_coords]).flatten()
    d = {"labels":np.asarray(label_list), "C_Dist":all_dist}
    distances = pd.DataFrame(data = d)
    ## change labels to ResID (three letter code and pose number) and seqID for checking at in pymol 
    all_labels = distances["labels"].values
    res_ids = [x.split(",")[0] for x in all_labels]
    distances["ResID"]=res_ids
    seq_ids = [x.split(",")[1] for x in all_labels]
    distances["seqID"] = seq_ids
    res_type = [x.split(",")[2] for x in all_labels]
    distances["type"]=res_type
    # remove anything that isn't a residue (i.e. metals)
    distances = distances.loc[distances['type']=="protein"]
    del distances['labels']

    #print(distances.head(10))
    return(distances) 


def get_rosetta_energy_terms(directory, site, this_protein, cutoff = 15):
    #print("\n\n\nWorking in features_EnergyTerms:")
    pdb_id=site.get_struc_id()
    site_id=site.metal_atoms[0].SITE_ID
    SITE_center_pt=site.get_center()

    by_res_scores = get_res_scores(directory, pdb_id)
    c_dist = calculate_dist_to_Cs(this_protein, SITE_center_pt)
    by_res_scores = pd.merge(by_res_scores, c_dist, on = "ResID")
    
    by_res_scores = by_res_scores.sort_values(['C_Dist'], ascending=[True])
    by_res_scores["C_Dist^2"] = by_res_scores["C_Dist"].astype(float)*by_res_scores["C_Dist"].astype(float)
    #print(by_res_scores[['ResID', 'seqID', 'resName', 'C_Dist', 'fa_atr', 'lk_ball', 'total']].iloc[:5])
    #print("using df with shape:", by_res_scores.shape)

    by_res_scores["tmp"] = 0
    energy_term_df=pd.DataFrame([site_id],["SITE_ID"]).T
    for term in energy_cols:
        by_res_scores["tmp"] = by_res_scores[term].astype(float)/by_res_scores["C_Dist^2"]
        energy_term_df["rbf_%s"%(term)] = by_res_scores.loc[by_res_scores['C_Dist']<=cutoff, "tmp"].sum()

    num_resis_used = len(by_res_scores.loc[by_res_scores['C_Dist']<=cutoff])
    energy_term_df["rbf_numScoredResis"] = num_resis_used
    del energy_term_df['SITE_ID']
    return(energy_term_df.T)

def get_rosetta_surface_terms(directory, site, this_protein, cutoff_list=[5,10,15,20,25,30, None]):
    pdb_id=site.get_struc_id()
    site_id=site.metal_atoms[0].SITE_ID
    SITE_center_pt=site.get_center()

    by_res_scores = get_res_scores(directory, pdb_id)
    c_dist = calculate_dist_to_Cs(this_protein, SITE_center_pt)
    by_res_scores = pd.merge(by_res_scores, c_dist, on = "ResID")
    by_res_scores = by_res_scores.sort_values(['C_Dist'], ascending=[True])

    surface_term_df=pd.DataFrame([site_id],["SITE_ID"]).T
    for term in ["BSA", "SASA"]:
        for cutoff in cutoff_list:
            if None!=cutoff:
                surface_term_df["SA%s_%s"%(cutoff, term)] = by_res_scores.loc[by_res_scores['C_Dist']<=cutoff, term].sum()
            else:
                surface_term_df["SA_All_%s"%(term)] = by_res_scores[term].sum()
    del surface_term_df['SITE_ID']
    return(surface_term_df.T)

