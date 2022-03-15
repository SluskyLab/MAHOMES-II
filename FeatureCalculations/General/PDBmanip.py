import numpy as np
import os
import glob
import decimal
import subprocess
import scipy.spatial
import math
import itertools

aa = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE", "SEC", "CH6"]
oneletAA = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "M", "C", "MYG"]
noncanAA = ["MSE", "IAS", "SEC"] #selenomethionene, L-asp acid for crosslinking, selenocysteine
metals = ["CUA", "CU", "FE", "MG", "ZN"]

def calc_RMSD(res_coords_native, res_coords_model):
    if len(res_coords_model) != len(res_coords_native):
        print("Uneven number of atoms", len(res_coords_native), len(res_coords_model))
        return
    rmsd = np.sqrt(np.sum((res_coords_native - res_coords_model)**2)/len(res_coords_native))
    return rmsd

def one_chain_pdb(filename, pdb_id, chainID = "A", keep_header = True, remove_tags = True, remove_RosLink = False, multiple_occupancy = True, inplace = True, outFile = "default"):
    #filename should be full filepath; pdb_id is 4-digit pdbID code, used for saving outFile
    #chainID should be a valid letter/number that exists in the pdb
    #keep_header == True will keep all lines before the ATOM lines
    #remove_tags skips residues that have been marked as expression tag, purification tag, initiating methionine, initiating residue or leader seuqence in the SEQADV lines
    #multiple_occupancy only keeps atoms that are coded with " ", "A", or "1" in the character position
    #inplace will write over the original file; if this equals false, the default file path will be the same directory as the original with PDB_Chain as the output file
    
    #identify residues with tags; 
    bad_res = []
    if remove_tags == True:
        bad_res = get_tag_resIDs(filename, chainID)
    #list of muliple occupancy codes to KEEP
    if multiple_occupancy == False:
        MO_codes = [" ", "A", 'B', 'C', 'D', 'E', '1', '2', '3', '4', '5']
    else:
        MO_codes = [" ", "A", "1"] 
    
    #determine where to save file
    if inplace == True:
        out_pdb = filename
    else:
        if outFile == "default": #if inplace == False, but no outFile has been designated, make one up from the pdb and chain
            file_path = filename.split("/")[:-1]
            #print(file_path)
            if len(file_path) > 0: 
                out_pdb = "/".join(file_path)
                out_pdb += "/%s_%s.pdb"%(pdb_id, chainID)
            else:
                out_pdb = "%s_%s.pdb"%(pdb_id, chainID)
        else: #or just use the one specified
            out_pdb = outFile
    
    #work through lines of file
    removed_res = False
    with open(filename, "r") as orig_pdb:
        orig_pdb = orig_pdb.readlines()        
    with open(out_pdb, "w+") as outData:
        atoms_reached = False
        for line in orig_pdb:
            if line[0:6] == "ENDMDL": #stop at the end of the first model for NMR structures
                outData.write(line)
                break
            
            if line[0:4] == "ATOM":
                #if line[16] != " ":
                    #print(line)
                atoms_reached = True
                if line[21] == chainID and line[22:26].strip() not in bad_res and line[16] in MO_codes:
                    #print(line[21], line)
                    outData.write(line)
                if line[22:26].strip() in bad_res:
                    removed_res = True
            elif line[0:6] == "HETATM" and line[21] == chainID and line[16] in MO_codes:
                outData.write(line)
            elif line[0:3] == "TER" and line[21] == chainID and line[16] in MO_codes:
                outData.write(line)
            elif atoms_reached == False and keep_header == True:
                if line[0:4] == "LINK" and remove_RosLink == True:
                    outData.write("#" + line)
                else:
                    outData.write(line)
    if removed_res == True:
        print("Removed residues from SEQADV lines of", pdb_id, chainID)

def one_chain_lines(orig_pdb_lines, chainID, multiple_occupancy = True):
    if multiple_occupancy == False:
        MO_codes = [" ", "A", 'B', 'C', 'D', 'E', '1', '2', '3', '4', '5']
    else:
        MO_codes = [" ", "A", "1"]
        
    keep_lines = []
    for line in orig_pdb_lines:
        if line[0:6] == "ENDMDL": #stop at the end of the first model for NMR structures
            break        
        if line[0:4] == "ATOM" and line[21] == chainID and line[16] in MO_codes:
            keep_lines.append(line)
        elif line[0:6] == "HETATM" and line[21] == chainID and line[16] in MO_codes:
            keep_lines.append(line)
        elif line[0:3] == "TER" and line[21] == chainID and line[16] in MO_codes:
            keep_lines.append(line)
    return keep_lines
    

def get_tag_resIDs(filename, chainID):
    tags = ["EXPRESSION TAG", "PURIFICATION TAG", "INITIATING METHIONINE", "INITIATING RESIDUE", "LEADER SEQUENCE"]
    try:
        seq_adv = subprocess.check_output(["grep", "^SEQADV", filename])
        seq_adv = seq_adv.decode("utf-8").strip().split("\n")
        #print(seq_adv)
        bad_res = []
        for line in seq_adv:
            #print(line[49:].strip(), line[16])
            if line[49:].strip() in tags and line[16] == chainID:
                bad_res.append(line[18:22].strip())
        #print(bad_res)
    except subprocess.CalledProcessError:
        #print("No SEQADV lines")
        bad_res = []
        
    return bad_res

def renumber_pdb_contig(filename, start_value = 1, inplace = True):
    with open(filename, "r") as orig_pdb:
        orig_pdb = orig_pdb.readlines()
    if inplace == True:
        new_filename = filename
    else: 
        new_filename = filename[:-4] + "Renumb.pdb"
    
    with open(new_filename, "w+") as new_pdb:    
        prev_res = -999
        count = start_value - 1
        for line in orig_pdb:
            if "ANISOU" in line: continue
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                this_res = int(line[22:26].strip())
                if this_res == prev_res:
                    new_res_num = (4-len(str(count)))*" " + str(count)
                    new_pdb.write(line[0:22] + new_res_num + line[26:])
                    prev_res = this_res
                else:
                    count += 1
                    new_res_num = (4-len(str(count)))*" " + str(count)
                    new_pdb.write(line[0:22] + new_res_num + line[26:])
                    prev_res = this_res
            else:
                new_pdb.write(line)
 
def renumber_pdb(pdb_file, filename, start_value):
    with open(filename, "w+") as new_pdb:
        count = 0
        adj_factor = 0
        for line in pdb_file:
            if "ANISOU" in line: continue
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                if count == 0:
                    count +=1
                    old_start = int(line[22:26].strip()) 
                    adj_factor = start_value - old_start
                    #print(adj_factor)
                    new_res = int(line[22:26].strip()) + adj_factor
                    new_res_num = (4-len(str(new_res)))*" " + str(new_res)
                    new_pdb.write(line[0:22] + new_res_num + line[26:])
                else:
                    new_res = int(line[22:26].strip()) + adj_factor
                    new_res_num = (4-len(str(new_res)))*" " + str(new_res)
                    new_pdb.write(line[0:22] + new_res_num + line[26:])
            else:
                new_pdb.write(line)
                
def renumber_pdb_inplace(pdb_file, start_value):
    new_lines = []
    count = 0
    adj_factor = 0
    for line in pdb_file:
        if "ANISOU" in line: continue
        if line[0:4] == "ATOM" or line[0:6] == "HETATM":
            if count == 0:
                #print(line)
                count +=1
                old_start = int(line[22:26].strip()) 
                adj_factor = start_value - old_start
                #print(adj_factor)
                new_res = int(line[22:26].strip()) + adj_factor
                new_res_num = (4-len(str(new_res)))*" " + str(new_res)
                new_lines.append(line[0:22] + new_res_num + line[26:])
            else:
                new_res = int(line[22:26].strip()) + adj_factor
                new_res_num = (4-len(str(new_res)))*" " + str(new_res)
                new_lines.append(line[0:22] + new_res_num + line[26:])
        else:
            new_lines.append(line)
    return(new_lines, adj_factor)
                                
def renumber_pdb_Rosetta(pdb_file, filename, start_value):
    #pdb_file is the readlines input
    with open(filename, "w+") as new_pdb:
        count = 0
        adj_factor = 0
        for line in pdb_file:
            if "ANISOU" in line: continue
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                if count == 0:
                    chain = line[21]
                    count +=1
                    old_start = int(line[22:26].strip()) 
                    adj_factor = start_value - old_start
                    #print(adj_factor)
                    new_res = int(line[22:26].strip()) + adj_factor
                    new_res_num = (4-len(str(new_res)))*" " + str(new_res)
                    new_pdb.write(line[0:22] + new_res_num + line[26:])
                elif line[21] == chain:
                    new_res = int(line[22:26].strip()) + adj_factor
                    new_res_num = (4-len(str(new_res)))*" " + str(new_res)
                    new_pdb.write(line[0:22] + new_res_num + line[26:])
                    end_res = new_res
                elif line[21] != chain:
                    chain = line[21]
                    old_start = int(line[22:26].strip()) 
                    adj_factor = end_res
                    #print(adj_factor)
                    new_res = int(line[22:26].strip()) + adj_factor
                    new_res_num = (4-len(str(new_res)))*" " + str(new_res)
                    new_pdb.write(line[0:22] + new_res_num + line[26:])
                    end_res = new_res
            else:
                new_pdb.write(line)

def ID_neighbors(res_list, resnum, target_res, cutoff):
    #this identifies neighbors within cutoff of ALL res in target_res list
    cutoff = int(cutoff)
    all_neighbors = []
    if (isinstance(target_res, list) == False): target_res = [target_res]
    for value in target_res:
        neighbors = []
        id_value = resnum.index(int(value))
        for res in res_list:
            if np.min(scipy.spatial.distance.cdist(res.Coords, res_list[id_value].Coords, "euclidean")) <= cutoff: 
                neighbors.append(res.resnum)
        if len(all_neighbors) == 0: 
            all_neighbors.extend(neighbors)
        else:
            all_neighbors = sorted(set(all_neighbors).intersection(neighbors))
        #print(all_neighbors)
    all_neighbors = list(sorted(set(all_neighbors).difference(target_res)))
    return(all_neighbors)
    
def ID_neighbors_perl(pdb1, target_res, cutoff):
    #this checks for any res with atom w/in cutoff of any target_res in list
    current_dir = os.getcwd()
    pdb_path = os.path.join(current_dir, pdb1)
    #print(type(target_res))
    if (isinstance(target_res, list) == True): target_res = ":".join(map(str, target_res))
    vicinity = subprocess.check_output(["/Users/meghan/mmtsb/perl/vicinity.pl", "-l", "%s" %target_res, "-hard", "%s" %cutoff, "-soft", "%s" %cutoff - 1, "%s.pdb" %pdb_path])
    vicinity = vicinity.decode("utf-8").strip()
    print(vicinity)
    if (len(vicinity) >2): chain = vicinity[0]
    else: chain = ""
    
    #return of vicinity.pl is in form of A32:45=A47=A52:55; "=" is separation between ranges
    vicinity = vicinity.split("=")    
    #print(chain)
    vicinity = [k[1:] for k in vicinity if k != "+1"]
    vicinity_list = []
    for value in vicinity:
        if ":" in value:
            values = value.split(":")
            x = 0
            for x in range(int(values[0]),int(values[1])+1): vicinity_list.append(str(x))
        else: vicinity_list.append(value)
    #print(vicinity_list)
    return(vicinity_list, chain)
            
def align_seq(seq1, seq2, pdb1 = "pdb1", pdb2 = "pdb2", clustalWpath = "/Users/meghan/CompBioPrograms/ClustalW2/clustalw2"):
    #there is a better version that is more customized in Sequences.py
    with open("%s_%sSeq.txt" %(pdb1, pdb2), "w+") as seq_file:
        seq_file.write(">" + pdb1 + "\n")
        seq_file.write(seq1 + "\n")
        seq_file.write(">" + pdb2 + "\n")
        seq_file.write(seq2 + "\n")
    clustalW = subprocess.check_output([clustalWpath, "-INFILE=%s_%sSeq.txt" %(pdb1, pdb2)])
    
    pdb1_al_seq = []
    pdb2_al_seq = []
    with open("%s_%sSeq.aln" %(pdb1, pdb2), "r") as aligned_seq:
        for line in aligned_seq:
            if pdb1 in line:
                line = line.strip().split()
                #print(line)
                pdb1_al_seq.extend(line[1])
            elif pdb2 in line:
                line = line.strip().split()
                #print(line)
                pdb2_al_seq.extend(line[1])
            else:
                continue
    #print(pdb1_al_seq)
    #print(pdb2_al_seq)
    os.system("rm *.aln")
    os.system("rm *.dnd")
    return pdb1_al_seq, pdb2_al_seq

def get_res_numbers(pdb1, chainID = "A"):
    res_list = []
    with open(pdb1, "r") as pdb_file:
        for line in pdb_file:
            if line[0:4] == "ATOM" and line[12:16].strip() == "CA"and line[21] == chainID:
                res_list.append(line[22:26].strip())
            elif line[0:6] == "HETATM" and line[12:16].strip() == "CA" and line[17:20] in aa and line[21] == chainID:
                res_list.append(line[22:26].strip())
            elif "ENDMDL" in line: break
    return res_list

def get_POR(pdb1, pdb2):
    pdb1_seq, pdb2_seq = align_seq(pdb1, pdb2)
    pdb1_res = get_res_numbers(pdb1)
    pdb2_res = get_res_numbers(pdb2)
    #print(pdb1_res)
    #print(pdb2_res)
    k = 0
    while k <= len(pdb1_seq)-1:
        if pdb1_seq[k] == "-":
            pdb1_res.insert(k, "-")
        if pdb2_seq[k] == "-":
            pdb2_res.insert(k, "-")
        k+= 1
    #print(pdb1_res)
    #print(pdb2_res)
    return pdb1_res, pdb2_res
    
def align_structures_MMTSB(pdb1, pdb1chain, pdb2):
    #this doesn't quite work...
    #mmtsb alligns by res number; therefore, residues must match
    #aligns SECOND pdb to first pdb; second pdb will be renumbered and its chain changed to match
    pdb1_seq, pdb2_seq = align_seq(pdb1, pdb2)
    pdb1_res = get_res_numbers(pdb1)
    pdb2_res = get_res_numbers(pdb2)
    #print(pdb1_res)
    #print(pdb2_res)
    k = 0
    while k <= len(pdb1_seq)-1:
        if pdb1_seq[k] == "-":
            pdb1_res.insert(k, "-")
        if pdb2_seq[k] == "-":
            pdb2_res.insert(k, "-")
        k+= 1
    #print(pdb1_res)
    #print(pdb2_res)
    with open("%s.pdb" %pdb2, "r") as pdb_file:
        old_pdb = pdb_file.readlines()
    with open("%s_renumb.pdb" %pdb2, "w+") as new_pdb:
        for line in old_pdb:
            if line[0:4] == "ATOM":
                res_num = line[22:26].strip()
                new_res = pdb1_res[pdb2_res.index(res_num)]
                new_pdb.write(line[0:21] + pdb1chain + new_res + line[26:])
    os.system("/Users/meghan/mmtsb/perl/lsqfit.pl -sel heavy %s.pdb %s_renumb.pdb 2>&1 1>%s_aligned.pdb" %(pdb1, pdb2, pdb2)) #-resnumonly
    #os.remove("%s_renumb.pdb" %pdb2)
    print("RMS of aligned structures is: ")
    os.system("/Users/meghan/mmtsb/perl/rms.pl %s.pdb %s_aligned.pdb" %(pdb1, pdb2))

def align_structures(pdb1, pdb2):
    #pdb2 will be aligned to pdb2
    os.system("/Users/meghan/TMalign -A %s.pdb -B %s.pdb >TM_%s.txt" %(pdb1, pdb2, pdb1))
    
    with open("TM_%s.txt" %pdb1, "r") as TM_align:
        tm_lines = TM_align.readlines()
    
    x_matrix = tm_lines[20].strip().split()
    x_matrix = x_matrix[1:]
    y_matrix = tm_lines[21].strip().split()
    y_matrix = y_matrix[1:]
    z_matrix = tm_lines[22].strip().split()
    z_matrix = z_matrix[1:]
    rot_matrix = np.zeros((3,3), dtype = float)
    k = 0
    for k in range(0, 3): 
        rot_matrix[0][k] = x_matrix[k+1]
        rot_matrix[1][k] = y_matrix[k+1]
        rot_matrix[2][k] = z_matrix[k+1]
    trans_matrix = [float(x_matrix[0]), float(y_matrix[0]), float(z_matrix[0])]
    #print(rot_matrix, trans_matrix)
    #rot_matrix[0][1] = x_matrix[1]
    
    with open("%s.pdb" %pdb2, "r") as pdb_file:
        old_pdb = pdb_file.readlines()
    
    X_coords = np.zeros(len(old_pdb), dtype = float)
    Y_coords = np.zeros(len(old_pdb), dtype = float)
    Z_coords = np.zeros(len(old_pdb), dtype = float)    
    for i in range(0, len(old_pdb)):
        if old_pdb[i][0:4] == "ATOM":
            X_coords[i] = float(old_pdb[i][30:38])
            Y_coords[i] = float(old_pdb[i][38:46])
            Z_coords[i] = float(old_pdb[i][46:54])
          
    all_coords = np.zeros((len(X_coords),3))
    j= 0
    for j in range(0, len(X_coords)):
        all_coords[j][0] = trans_matrix[0] + rot_matrix[0][0]*X_coords[j] + rot_matrix[0][1]*Y_coords[j] + rot_matrix[0][2]*Z_coords[j]
        all_coords[j][1] = trans_matrix[1] + rot_matrix[1][0]*X_coords[j] + rot_matrix[1][1]*Y_coords[j] + rot_matrix[1][2]*Z_coords[j]
        all_coords[j][2] = trans_matrix[2] + rot_matrix[2][0]*X_coords[j] + rot_matrix[2][1]*Y_coords[j] + rot_matrix[2][2]*Z_coords[j]
    
    #print("%0.3f" %all_coords[0][0]) # decimal.Decimal.from_float(all_coords[0][0]))
    with open("%s_aligned.pdb" %pdb2, "w+") as new_pdb:
        for i in range(0, len(old_pdb)):
            line = old_pdb[i]
            x_coord = "%0.3f" % all_coords[i][0]
            y_coord = "%0.3f" % all_coords[i][1]
            z_coord = "%0.3f" % all_coords[i][2]
            x_coord = x_coord.rjust(8)
            y_coord = y_coord.rjust(8)
            z_coord = z_coord.rjust(8)
            coords = x_coord+y_coord+z_coord
            #print(len(coords))
            new_pdb.write(line[0:30] + coords + line[54:])
    
    