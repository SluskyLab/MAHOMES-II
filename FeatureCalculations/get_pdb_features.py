import sys
import pandas as pd

import features_geom as f_geom
import features_EnergyTerms as f_energy
import features_electro as f_electro
import features_pocket as f_pocket
import Reformat

RESOURCE_DIR = "General/"
sys.path.insert(0, "%s" % RESOURCE_DIR)
import SITE as SITE
import PDBparser as pdbp

metal_size = {
         'MO':1, 'MOO':5, '4MO' : 1,'6MO' :1,'MOS': 4,
         'MG': 1,'3NI':1,'NI' : 1, 'ZN': 1,'MGF': 4,'MN3' : 1,'MN' : 1,'CO': 1,
         'OFE': 2, 'FE2': 1,'FEO': 3, 'FE' : 1,'FES': 4,
         'CU': 1, 'C2O' :3, 'CUA' : 2, 'CU1': 1, "3CO": 1,
         }


def get_site_info(struc_id, site):
    # different ways to count a site for ion codes, number of metals, and number of any atom (i.e. FES has 1, 2, and 4)
    site_atom_count = 0
    unique_seqNums= []; cur_i=0
    for metal_atom in site.metal_atoms:
        if not (metal_atom.seqID in unique_seqNums):
            unique_seqNums.append(metal_atom.seqID)
            if metal_atom.resName in metal_size.keys():
                site_atom_count += metal_size[metal_atom.resName]
        cur_i+=1

    SITE_info = pd.DataFrame([struc_id, site.SITE_ID, len(unique_seqNums), len(site.metal_atoms), site_atom_count], index=["struc_id", "SITE_ID", "MetalCodes", "MetalAtoms", "SiteAtoms"])
    SITE_info.loc["bad_site"] = site.bad
    SITE_info.loc["note"] = site.note

    cur_metal = 1
    ## add site metal atom info
    for metal in site.metal_atoms:
        SITE_info.loc["metal%s_atom_name"%str(cur_metal)] = metal.atom_name 
        SITE_info.loc["metal%s_serial"%str(cur_metal)] = metal.serial 
        SITE_info.loc["metal%s_resName"%str(cur_metal)] = metal.resName 
        SITE_info.loc["metal%s_seqID"%str(cur_metal)] = metal.seqID 
        cur_metal+=1
    return(SITE_info)

def get_features_for_site(struc_id, struc_dir, site, this_protein, res_nums, old_features=False):
    SITE_info=get_site_info(struc_id, site)
    ######################################
    #####    coordination geometry   #####  
    ######################################
    ## make identifiers for findgeo results
    findgeo_metals = []
    for metal in site.metal_atoms:
        new_metal = "%s_%s_%s_%s"%(metal.atom_name, metal.seqID, metal.serial, metal.chainID)
        findgeo_metals.append(new_metal)
    #print(findgeo_metals)
    geom_feat = f_geom.collate_findgeo(findgeo_metals, struc_id, struc_dir)
    #print(geom_feat)

    ######################################
    #####   Rosetta Energy Terms     #####  
    ######################################
    energy_terms = f_energy.get_rosetta_energy_terms(struc_dir, site, this_protein)
    #print(energy_terms)

    ######################################
    #####       electrostatics       #####  
    ######################################
    try:
        bluues_electro = f_electro.get_electro_features(struc_dir, site, this_protein)
    except:
        bluues_electro=[]
        site.bad = True
        site.note = "Failed bluues_electro feature calculations"
        #print("\tFailed bluues_electro %s site with %s %s"%(struc_id, site.metal_atoms[0].resName, str(site.metal_atoms[0].seqID)))
        SITE_info.loc["bad_site"] = site.bad
        SITE_info.loc["note"] = site.note
    #print(bluues_electro)

    ######################################
    #####        ghecom pocket       #####  
    ######################################
    pocket_shape, pocket_lining, pocket_info = f_pocket.pocket_lining_features(struc_dir, site, this_protein, res_nums)
    #print(pocket_shape)
    
    ######################################
    ####  combine feature domains    #####  
    ######################################
    calc_feats = []
    for domain in [SITE_info, bluues_electro, geom_feat, energy_terms, pocket_shape, pocket_lining, pocket_info]:
        #print(domain, len(domain))
        if len(domain)>0:
            calc_feats.append(domain)
    
    #dump the whole thing into a single dataframe and append to growing list
    site_features = pd.concat(calc_feats).T
    return(site_features)


def get_features_for_pdb(job_dir, struc_id, old_features=False):
    struc_features = []
    struc_dir = "%s/%s"%(job_dir, struc_id)
    struc_file = "%s/%s.pdb"%(struc_dir, struc_id)

    try:
        Reformat.Rosetta(struc_dir, struc_id)
        residues, res_nums, header = pdbp.create_res(struc_file)
        this_protein = pdbp.Protein(residues, res_nums, header)
        sites = SITE.find_SITEs(struc_id, struc_file)
    except:
        print("No PDB file output from Rosetta", struc_id)
        return_df = pd.DataFrame([[struc_id, "unknown", True, "Invalid structure file input. Please make sure file is valid input for Rosetta 3.13."]], columns=["struc_id", "SITE_ID", "bad_site", "note"])
        return(return_df)

    bluues_output=True
    try:
        Reformat.bluues(struc_dir, struc_id)
    except:
        bluues_output=False
        print("No bluues output", struc_id)

    for site in sites:
        try:
            if True==site.bad:
                site_features =get_site_info(struc_id, site).T
            elif False==bluues_output:
                site.bad = True
                site.note = "Missing bluues"
                site_features =get_site_info(struc_id, site).T
            else:
                site_features = get_features_for_site(struc_id, struc_dir, site, this_protein, res_nums, old_features)
        except:
            site.bad = True
            site.note = "Failed feature calculations"
            print("\tFailed %s site with %s %s"%(struc_id, site.metal_atoms[0].resName, str(site.metal_atoms[0].seqID)))
            site_features = get_site_info(struc_id, site).T
                
        site_features.set_index("SITE_ID", drop=True, inplace=True)
        struc_features.append(site_features)

    if 1<len(struc_features):
        struc_features = pd.concat(struc_features, join='outer', axis=0)
    elif 1==len(struc_features):
        struc_features=struc_features[0]
    elif 0==len(struc_features):
        return(pd.DataFrame())
    struc_features.reset_index(drop=False, inplace=True)
    return(struc_features)



def run_sample(old_features=True):
    #test_dir = "/Users/ryan/Desktop/data/MAHOMESII/sample_data"
    #test_dir = "/Users/ryan/Desktop/data/MAHOMESII/dataset_calculations/data"
    #test_pdb ="1ard_A"

    test_dir = "/Users/ryan/Desktop/data/MAHOMESII/mutant_calculations/data"
    test_pdb ="4k2s_H_Relax1"

    test_features = get_features_for_pdb(test_dir, test_pdb, old_features)
    print(test_features.shape)
    print(test_features)
    features = test_features.columns.tolist()
    features.sort()
    for term in features:
        print(term, test_features[term].iloc[0])
#run_sample(old_features=False)
