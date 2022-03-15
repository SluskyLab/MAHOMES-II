import numpy as np
import os
import glob
import decimal
import subprocess
import scipy.spatial
import math
import itertools
import PDBmanip as pdbm

aa = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE", "SEC"]
oneletAA = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "M", "C"]
mod_res_list = ["ABA", "DDZ", #modified residues are generally noncanonical and/or chemically modified for linkages, etc but they are much less reconizable than simple swaps of S -> Se
                "SSN", 
                "BFD", "PHD", "IAS",
                "CAF", "CAS", "CME", "CMH", "CSD", "CSO", "CSS", "CSX", "OCS", "QCS", "SMC", "SNC", "YCM", 
                "PCA",
                "HIC", "HIP", "HS8", "NEP", "OHI",
                "6CL", "ALY", "BTK", "FAK", "KCR", "KCX", "LLP", "M3L", "MLZ", "SLL",
                "FME", "MHO", "MME", "MSO",
                "PHA",
                "05N", "HYP", "PXU",
                "SEP", 
                "TPO",
                "TRX", "NIY",
                "OTY", "PTR", "TY2",
                "MVA"]
mod_res_aa = ["A", "A", "N", "D", "D", "D", "C", "C", "C", "C", "C", "C", "C", "C", "C","C","C", "C","C", #this is the closest parent amino acid to the modified residue
            "Q", "H", "H", "H", "H", "H", "K", "K", "K", "K", "K", "K", "K", "K", "K", "K", 
            "M", "M", "M", "M", "F", "P", "P", "P", "S", "T", "W", "W", "Y", "Y", "Y", "V"]
chromophore_res = ["DYG", "CRQ", "CRG", "CRO", "GYS", "OFM", "CRG"]
metals = ["CUA", "CU", "FE", "MG", "ZN", "MN", "FE2", "FES", "FEO", "CU1", "MO", "MOO", "MOS", "NI", "3CO", "CO", "K", "CA", "NA"]
header_delims = ["HEADER", "SEQRES", "HET   ", "HETNAM", "EXPDTA", "SOURCE", "COMPND", "TITLE ", "SEQADV", "MODRES"] #this should include other header start info of relevance with 6 characters

#a Protein is created by the following pair of calls:

#residues, res_nums, header = pdbp.create_res("1yew.pdb")
#this_protein = pdbp.Protein(residues, res_nums, header)

#this produces an object Protein with numerous properties of use to our work
#the list Protein.residues are class Residues with their own properties; these include all atoms and hetatom lines

class Protein:
    def __init__(self, res_list, res_nums, header_info):
        #print(header_info)
        self.sequence = "" #from atom record lines
        self.NumRes = len(res_list)
        self.num_atoms = 0
        self.Coords = []
        self.Coords_index = []
        self.metals = []
        self.residues = res_list
        self.res_nums = []
        self.gene_seq = {} #from SEQRES lines by chain
        self.chains = [] #sorted set of chains associated with residues
        self.mutated = False
        self.chromophore = False
        self.modres = [] #res numbers of modified residues
        self.waters = [] #res numbers of waters; water does not get added to the coord list - don't want it for neighbor-finding purposes
        
        for x in res_list:
            self.sequence += x.onelet
            self.num_atoms += len(x.Atoms)
            new_resnum = str(x.resnum) + x.chain
            self.chains.append(x.chain)
            self.res_nums.append(new_resnum)
            if x.type == "metal": self.metals.append(new_resnum)
            elif x.type == "water": self.waters.append(new_resnum)
            elif x.modres == True: self.modres.append(new_resnum)
            elif x.type == "chromophore": self.chromophore = True
            
            if x.type != "water":
                self.Coords.extend(x.Coords)
                for atom in x.Atoms:
                    self.Coords_index.append(new_resnum)

        self.chains = sorted(set(self.chains))
        for x in self.chains:
            self.gene_seq[x] = ""
        self.Coords = np.asarray(self.Coords)
        self.KDTree = scipy.spatial.KDTree(self.Coords)
        
        for line in header_info:
            if "SEQRES" in line:
                line = line.split()
                #print(line)
                for value in line[4:]:    
                    try:
                        #print(value)
                        self.gene_seq[line[2]] += oneletAA[aa.index(value)]
                    except ValueError:
                        self.gene_seq[line[2]] += "X"
                    except KeyError:
                        continue
            if "MUTATION: YES" in line:
                self.mutated = True
        
    def get_neighbor_res(self, atom_coords, r):
        neighbor_atoms = self.KDTree.query_ball_point(atom_coords, r)
        if np.ndim(atom_coords) > 1: 
            neighbor_atoms = [y for x in neighbor_atoms for y in x]
        #for x in neighbor_atoms:
            #print(self.Coords[x], self.Coords_index[x], np.sqrt(sum(np.square(atom_coords - self.Coords[x]))) )
        neighbor_res = [self.Coords_index[x] for x in neighbor_atoms]
        neighbor_res = sorted(set(neighbor_res))
        return neighbor_res
        
    def get_metal_shells(self, metals, cutoff1, cutoff2):
        print("There are this many metals:", len(metals), metals)
        metal_first_shell = [[] for x in metals]
        for x in range(0, len(metals)):
            first_shell = self.get_neighbor_res(self.residues[self.res_nums.index(metals[x])].Coords, cutoff1)
            #print(metals[x], first_shell)
            metal_first_shell[x] = first_shell #sorted(set(first_shell).difference([metals[x]])) #don't actually want to remove metal for the mononuclear cases
            #print(metals[x], metal_first_shell[x])
        #print(metal_first_shell)
        
        metal_sites = []
        for i, x in enumerate(metal_first_shell):
            #print("\n", i, x)
            for y in range(len(metals)):
                #check if metals are w/in 3A, then check if any res are in common between metals
                if metals[y] in x: 
                    metal_first_shell[y].extend(x)
                    metal_first_shell[i].extend(metal_first_shell[y])
                elif len(set(x).intersection(metal_first_shell[y])) != 0:
                    metal_first_shell[y].extend(x)
                    metal_first_shell[i].extend(metal_first_shell[y])
            metal_first_shell = [sorted(set(x)) for x in metal_first_shell]
            #print(metal_first_shell)
        metal_first_shell.sort()
        #print("finished first shell")
        #print(metal_first_shell)
        metal_first_shell = list(metal_first_shell for metal_first_shell, _ in itertools.groupby(metal_first_shell))
        for value in metal_first_shell:
            #print(value)
            metal_sites.append([sorted(set(metals).intersection(value)) , value])
        
        #print("Metal ligands are: ", metal_sites)
        for site in metal_sites:
            metal = site[0]
            site_res = site[1]
            #print(metal, site_res)
            site_coords = [self.residues[self.res_nums.index(x)].Coords for x in site_res]
            site_coords = [y for x in site_coords for y in x]
            site_coords = np.asarray(site_coords)    
            new_res = self.get_neighbor_res(site_coords, cutoff2)
            site[1] = new_res
            #pymol_res = sorted([int(x[:-1]) for x in new_res])
            #print("+".join(map(str, pymol_res)))
        return metal_sites #,pymol_res
                    
class Residue:
    def __init__(self, resnum, name, chain, atom_list, res_index):
        #print("atomlist is ", len(atom_list))
        self.name = name
        self.resnum = resnum #from PDB
        self.res_index = res_index #total count
        self.modres = False
        self.bfactors = []
        self.occupancy = []
        self.onelet = ""
        if self.name in metals:
            self.type = "metal"
        elif self.name == "HOH":
            self.type = "water"
        elif self.name in chromophore_res:
            self.type = "chromophore"
        elif self.name in aa: #modified 4/24/19 because crystal structures often contain SO4, BME, etc that aren't really residues and should be labeled other
            self.type = "protein"
            self.onelet = oneletAA[aa.index(self.name)]
        elif self.name in mod_res_list:
            self.type = "protein"
            self.modres = True
            self.onelet = mod_res_aa[mod_res_list.index(self.name)]
        else:
            self.type = "cofactor"
            
        self.chain = chain
        self.Atoms = []
        self.Coords = np.zeros((len(atom_list), 3))
        self.num_atoms = len(atom_list)
        self.psi = 500.00
        self.phi = 500.00
        self.chi1 = 500.00   
        for x in range(0, len(atom_list)):
            #print(atom_list[x])
            self.Atoms.append(atom_list[x][0])
            self.Coords[x][0] = atom_list[x][1]
            self.Coords[x][1] = atom_list[x][2]
            self.Coords[x][2] = atom_list[x][3]
            self.occupancy.append(atom_list[x][4])
            self.bfactors.append(atom_list[x][5])

    def set_phi(self, prev_res):
        #phi is prevC-N-C-C
        if prev_res.resnum + 1 == self.resnum and self.name not in metals:
            self.phi = dihedral(prev_res.Coords[2], self.Coords[0], self.Coords[1], self.Coords[2])
        
    def set_psi(self, next_res):
        #psi is N-c-c-nextN
        if self.resnum + 1 == next_res.resnum and self.name not in metals:
            self.psi = dihedral(self.Coords[0], self.Coords[1], self.Coords[2], next_res.Coords[0])
        
def chi1_4th(res_id):
    if (res_id == "VAL" or res_id == "ILE"): return "CG1"
    if res_id == "CYS": return "SG"
    if res_id == "SER": return "OG"
    if res_id == "THR": return "OG1"
    else: return "CG"
        
def dihedral(a, b, c, d):
    p = np.vstack((a, b, c, d))
    #print(p)
    # Calculate vectors between points, b1, b2, and b3 in the diagram
    b = p[:-1] - p[1:]
    #print(b)
    # "Flip" the first vector so that eclipsing vectors have dihedral=0
    b[0] *= -1
    # Use dot product to find the components of b1 and b3 that are not
    # perpendicular to b2. Subtract those components. The resulting vectors
    # lie in parallel planes.
    v = np.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
    # Normalize vectors
    v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
    #print(v)
    b1 = b[1] / np.linalg.norm(b[1])
    #print(b1)
    x = np.dot(v[0], v[1])
    #print(x)
    m = np.cross(v[0], b1)
    y = np.dot(m, v[1])
    return np.degrees(np.arctan2( y, x ))

def set_all_phi_psi(residues, res_nums):
    for value in range(0,len(res_nums)):
        if value == 0:
            residues[value].set_psi(residues[value+1])
        elif value == len(res_nums) -1:
            residues[value].set_phi(residues[value-1])
        else:
            #print("prev ", residues[value-1].Coords)
            #print("next ", residues[value+1].Coords)
            residues[value].set_phi(residues[value-1])
            residues[value].set_psi(residues[value+1])

class Atom:
    def __init__(self, name, resnum, x, y, z):
        self.name = name
        self.resnum = int(resnum)
        self.coords = [float(x), float(y), float(z)]

def convert_lines_to_res(pdb_lines):
    all_res = []
    res_nums = []
    header = []
    prev_res = -1000
    res_atoms = []
    resnum = 0
    chain = ""
    name = ""
    res_index = 0
    for line in pdb_lines:
        if line[0:4] == "ATOM" or line[0:6] == "HETATM":
            #print(res_index)
            new_res = int(line[22:26])
            #print(prev_res, resnum, len(res_atoms))
            if prev_res == new_res:
                occ = line[54:60].strip()
                try: occ = float(occ) 
                except ValueError: occ = np.nan
                bfact = line[60:66].strip()
                try: bfact = float(bfact) 
                except ValueError: bfact = np.nan
                #res_atoms.append([line[12:16].strip(), float(line[30:38]), float(line[38:46]), float(line[46:54]), float(line[54:60].strip()), float(line[60:66].strip())])
                res_atoms.append([line[12:16].strip(), float(line[30:38]), float(line[38:46]), float(line[46:54]), occ, bfact])
                #print(res_atoms)
            else:
                #print("atoms to append: ", len(res_atoms), name)
                if len(res_atoms) > 0:
                    #print("atoms present")
                    #print(res_atoms)
                    all_res.append(Residue(resnum, name, chain, res_atoms, res_index))
                    #print(resnum, name, chain, res_index)
                    res_nums.append(resnum)
                    res_index += 1
                res_atoms = []
                prev_res = int(line[22:26])
                resnum = int(line[22:26])
                name = line[17:20].strip()
                chain = line[21]
                occ = line[54:60].strip()
                try: occ = float(occ) 
                except ValueError: occ = np.nan
                bfact = line[60:66].strip()
                try: bfact = float(bfact) 
                except ValueError: bfact = np.nan
                #res_atoms.append([line[12:16].strip(), float(line[30:38]), float(line[38:46]), float(line[46:54]), float(line[54:60].strip()), float(line[60:66].strip())])
                res_atoms.append([line[12:16].strip(), float(line[30:38]), float(line[38:46]), float(line[46:54]), occ, bfact])
        elif line[0:6] in header_delims:
            header.append(line) 
    if len(res_atoms) != 0:
        #print(res_atoms) 
        all_res.append(Residue(resnum, name, chain, res_atoms, res_index))
        res_nums.append(resnum)
    return(all_res, res_nums, header)      
            
def create_res(pdb, lines_only = False):
    if lines_only == False:
        with open(pdb, "r") as pdb_file:
            pdb_file = pdb_file.readlines()
        these_res, these_nums, this_header = convert_lines_to_res(pdb_file)
    else:
        these_res, these_nums, this_header = convert_lines_to_res(pdb)
    return(these_res, these_nums, this_header)

#functions to retrieve PDB files
def download_pdbs(pdb_list, output_path = "", header_only = False):
    #pdb_list should be a list of 4-digit PDB_IDs; if not, it gets converted to a list first
    #if pdb_IDs aren't 4-digits: downloads pdb corresponding to first 4 digits - useful if you later need to extract a single chain
    #output_path should include the final /; the default is the current directory
    #header files can be retrieved with the optional flag header_only set to True
    if isinstance(pdb_list, list) == False:
        pdb_list = [pdb_list]
    for entry in pdb_list:
        print(entry)
        if header_only == False:
            #downloads PDB files and checks briefly; if there are less than 100 lines, the pdb is removed and the cif file retrieved instead
            if os.path.isfile("%s%s.pdb"%(output_path, entry[0:4])) == False:
                download_file(entry[0:4], output_path, "https://files.rcsb.org/download/")
            else:
                print("Already downloaded")
        if header_only == True:
            if os.path.isfile("%s%s.pdb"%(output_path, entry[0:4])) == False:
                download_file(entry[0:4], output_path, "https://files.rcsb.org/header/")
            else:
                print("Already downloaded header")

def download_file(pdb, output_path, database_path):
    if os.path.isfile("%s%s.cif"%(output_path, pdb)) == True:
        print("Already downloaded CIF file")
    else:
        os.system("curl %s/%s.pdb -o %s%s.pdb" %(database_path, pdb, output_path, pdb))
        with open("%s%s.pdb"%(output_path, pdb), "r") as inData:
            inData = inData.readlines()
        if len(inData) < 20:
            print("This is too big to be in pdb format!!")
            os.system("rm %s%s.pdb"%(output_path, pdb))
            if os.path.isfile("%s%s.cif"%(output_path, pdb)) == False:
                os.system( "curl %s%s.cif -o %s%s.cif" %(database_path, pdb, output_path, pdb) )

#functions to clean up PDBs; esp useful for working with Rosetta < version 3.7    
def remove_H(pdb1):
    with open("%s.pdb" %pdb1, "r") as origPDB:
        pdb_lines = origPDB.readlines()
    
    with open("%s.pdb" %pdb1, "w+") as new_PDB:
        for line in pdb_lines:
            if line[76:78].strip() != "H": new_PDB.write(line)

def remove_het(pdb1):
     with open("%s.pdb" %pdb1, "r") as origPDB:
         pdb_lines = origPDB.readlines()
    
     with open("%s.pdb" %pdb1, "w+") as new_PDB:
         for line in pdb_lines:
             if line[0:6] != "HETATM": new_PDB.write(line)
             
def clean_pdb(pdb1, HET = False, NoH = True):
    #This will remove any header info as well as any HET atoms by default. If Het is set to True, any HETATMs will be included in clean PDB
    with open("%s.pdb" %pdb1, "r") as origPDB:
        pdb_lines = origPDB.readlines()
    print(pdb_lines[0])
    with open("%s_clean.pdb" %pdb1, "w+") as newPDB:
        finished = False
        for line in pdb_lines:
            if line[0:6] == "EXPDTA":
                if "NMR" in line:
                    print("NMR structure")
                    finished = clean_NMR(pdb_lines, pdb1, HET, NoH)
                    #print(finished)
                else:
                    print("Not NMR")
                    finished = clean_xray(pdb_lines, pdb1, HET, NoH)
            elif finished == True:
                continue
            else: continue
            
def clean_xray(pdb_file, pdb1, het, NoH):
    with open("%s_clean.pdb" %pdb1, "w+") as newPDB:
        for line in pdb_file:
            if het == False:
                if line[0:4] == "ATOM":
                    if NoH == True and line[76:78].strip() != "H":
                        newPDB.write(line)
                    else:
                        newPDB.write(line)
            else:
                if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                    if NoH == True and line[76:78].strip() != "H":
                        newPDB.write(line)
                    else:
                        newPDB.write(line)
    return True

def clean_NMR(pdb_file, pdb1, het, NoH):
    #Clean NMR will only keep the first model of 20
    model = True
    with open("%s_clean.pdb" %pdb1, "w+") as newPDB:
        for line in pdb_file:
            #print(line[0:6], model)
            if line[0:6] == "ENDMDL":
                model = False
            elif model == True:
                if het == False and line[0:4] == "ATOM":
                    if NoH == True and line[76:78].strip() != "H":
                        newPDB.write(line)
                    else:
                        newPDB.write(line)
                elif het == True:
                    if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                        if NoH == True and line[76:78].strip() != "H":
                            newPDB.write(line)
                        else:
                            newPDB.write(line)
                else: continue
            else:
                continue
    return True
    
#miscellaneous, probably outdated functions relating to coordinates of PDBs; these are not written to be compatible with the class Protein but they probably should be
def get_CA_coords(pdb, res_list = [], chainID = "A"):
    with open("%s.pdb" %pdb, "r") as coord_file:
        if len(res_list) > 0:
            res_coords = np.zeros((len(res_list), 3))
            for line in coord_file:
                if line[12:16].strip() == "CA" and int(line[22:26]) in res_list and line[21] == chainID:
                    res_index = res_list.index(int(line[22:26]))
                    res_coords[res_index][0] = float(line[30:38])
                    res_coords[res_index][1] = float(line[38:46])
                    res_coords[res_index][2] = float(line[46:54])
        else:
            res_coords = []
            count = 0
            for line in coord_file:
                if line[12:16].strip() == "CA" and line[21] == chainID:
                    res_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    count +=1
                    #res_coords.extend(coord_line)
            res_coords = np.asarray(res_coords)
            #print(res_coords)
            res_coords = np.reshape(res_coords, (count, 3))
    
    return res_coords
    
def get_all_coords(pdb, res_list = [], chainID = "A"):
    res_list = [int(x) for x in res_list]
    with open("%s.pdb" %pdb, "r") as coord_file:
        if len(res_list) > 0:
            res_coords = np.zeros((len(res_list), 3))
            for line in coord_file:
                if int(line[22:26]) in res_list and line[21] == chainID:
                    res_index = res_list.index(int(line[22:26]))
                    res_coords[res_index][0] = float(line[30:38])
                    res_coords[res_index][1] = float(line[38:46])
                    res_coords[res_index][2] = float(line[46:54])
        else:
            res_coords = []
            count = 0
            for line in coord_file:
                if line[21] == chainID:
                    res_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    count +=1
                    #res_coords.extend(coord_line)
            res_coords = np.asarray(res_coords)
            print(res_coords)
            res_coords = np.reshape(res_coords, (count, 3))

    return res_coords

def get_res_centroids(pdb1, res_list = [], chainID = "A"):
    res_list = [int(x) for x in res_list]
    with open("%s.pdb" %pdb1, "r") as coord_file:
        if len(res_list) > 0:
            res_centroids = np.zeros((len(res_list),3))
            res_coords = []
            count = 0
            res_num = res_list[0]
            for line in coord_file:
                if line[0:4] == "ATOM" and line[12:16].strip() == "N" and len(res_coords) != 0:
                    res_coords = np.asarray(res_coords)
                    #print(res_coords)
                    res_coords = np.reshape(res_coords, (count, 3))
                    res_centroids[res_list.index(res_num)] = np.average(res_coords, axis = 0)
                    #print(res_coords)
                    #print(res_centroids[res_list.index(res_num)])
                    count = 0
                    res_coords = []
                
                res_num = int(line[22:26].strip())    
                if res_num in res_list and line[21] == chainID:
                    count += 1
                    res_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        else:
            res_centroids = []
            res_count = 0
            res_coords = []
            count = 0
            res_num = res_list[0]
            for line in coord_file:
                if line[0:4] == "ATOM" and line[12:16].strip() == "N" and len(res_coords) != 0:
                    res_coords = np.asarray(res_coords)
                    #print(res_coords)
                    res_coords = np.reshape(res_coords, (count, 3))
                    res_centroids.append(np.average(res_coords, axis = 0))
                    count = 0
                    res_count += 1
                    res_coords = []
                
                res_num = int(line[22:26].strip())    
                if line[21] == chainID:
                    count += 1
                    res_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    
            res_centroids = np.asarray(res_centroids)
            #print(res_centroids)
            res_centroids = np.reshape(res_centroids, (res_count, 3))

    return res_centroids
    