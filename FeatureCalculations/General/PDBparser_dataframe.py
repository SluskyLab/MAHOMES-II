#########################################################################################################################
#################################################    who (P)DB          #################################################
#########################################################################################################################

## file with general functions and variables useful for a lot of MSEAL stuff
import pandas as pd
import numpy as np
#import os
pd.set_option('display.max_colwidth', 100)
pd.set_option('display.max_rows', 10000)

#####################################################################
######             CREATE WORKABLE PDB DATAFRAME               ######
#####################################################################


## simple dataframe from pdb (with NMR adjustment)
def get_PDB_DF(pdb_id, PDB_file):
    f = open(PDB_file, "r")
    pdb_lines = f.readlines()
    pdb_file_df = pd.DataFrame(pdb_lines, columns = ['raw'])
    pdb_file_df['raw'] = pdb_file_df['raw'].str[:80]
    pdb_file_df['key'] = pdb_file_df['raw'].str[:6].str.replace(' ', '')
    ## only keep first model, mostly for NMR models
    pdb_file_df = get_first_model(pdb_file_df)
    
    return(pdb_file_df)


def get_first_model(PDB_DF):
    MODEL_DF = PDB_DF[PDB_DF['key'] == 'MODEL']
    ## return if NMR does not use multiple models
    if len(MODEL_DF) < 2:
        return(PDB_DF)
    
    ENDMDL_DF = PDB_DF[PDB_DF['key'] == 'ENDMDL']
    ## get range of all non first model coordinate lines
    index_srt_scnd_MDL = MODEL_DF.index.tolist()[1]
    index_end_MDLs = ENDMDL_DF.index.tolist()[len(ENDMDL_DF)-1]
    ## remove all non first model coordinate lines
    newPDB_DF = PDB_DF.loc[:index_srt_scnd_MDL-1]
    newPDB_DF = newPDB_DF.append(PDB_DF.loc[index_end_MDLs+1:])
    newPDB_DF = newPDB_DF.reset_index(drop=True)
    return(newPDB_DF)


#####################################################################
######  FOR VARIABLES THAT APPLY ON THE RESIDUE/ATOM LEVEL     ######
#####################################################################
## for help on any future work, use -URL->
## http://www.bmsc.washington.edu/CrystaLinks/man/pdb/guide2.2_frame.html

## 'ATOM' and 'HETATM' use get_ATOM_DF for a atom level dataframe
def get_ATOM_DF(ATOM, pdb=None):
    atom_data = split_ATOM(ATOM['raw'])
    ATOM = pd.DataFrame({
         'serial':atom_data[0]  , 'atom_name':atom_data[1]      ,  'altLoc':atom_data[2]     \
        ,'resName':atom_data[3] , 'chainID':atom_data[4]   ,  'seqNum':atom_data[5]     \
        ,'iCode':atom_data[6]   , 'x':atom_data[7]         ,  'y':atom_data[8]          \
        ,'z':atom_data[9]       , 'occupancy':atom_data[10] ,  'tempFactor':atom_data[11] \
        ,'segID':atom_data[12]   , 'element':atom_data[13]   ,  'charge':atom_data[14]     })
    
    #ATOM.loc[:,'Res_ID'] = ATOM['resName'].str[:] + ATOM['seqNum'].astype(str) + "_" + ATOM['chainID'].str[:]
    
    if pdb!=None:
        ## add mapping variables
        ATOM['RES_ID'] = pdb[:4] + "_" + ATOM['chainID'].str[:] + "_" +  ATOM['resName'].str[:] + "_" + ATOM['seqNum'].str[:]
        ATOM['ATOM_ID'] =  ATOM['RES_ID'].str[:] + "_" + ATOM['element'].str[:]+ "_" + ATOM['serial'].str[:] 
    # change coordinates to float to make them usable
    ATOM['x']=ATOM['x'].astype(float);ATOM['y']=ATOM['y'].astype(float);ATOM['z']=ATOM['z'].astype(float)
    return(ATOM)


def split_ATOM(raw):
    return(  raw.str[6:11].str.replace(' ', ''),  raw.str[12:15].str.replace(' ', ''),     raw.str[16].str.replace(' ', ''), \
            raw.str[17:20].str.replace(' ', ''),     raw.str[21].str.replace(' ', ''),  raw.str[22:26].str.replace(' ', ''), \
               raw.str[27].str.replace(' ', ''),  raw.str[30:37].str.replace(' ', ''),  raw.str[38:45].str.replace(' ', ''), \
            raw.str[46:53].str.replace(' ', ''),  raw.str[54:59].str.replace(' ', ''),  raw.str[60:65].str.replace(' ', ''), \
            raw.str[72:75].str.replace(' ', ''),  raw.str[76:78].str.replace(' ', ''),  raw.str[79:].str.replace(' ', ''))


## 'MODRES' lines describe how residues have been altered
def get_MODRES_DF(MODRES):
    modres_data = split_MODRES(MODRES['raw'])
    MODRES = pd.DataFrame({ 'idCode':modres_data[0], 'resName':modres_data[1], 'chainID':modres_data[2]  \
                           , 'seqNum':modres_data[3], 'iCode':modres_data[4], 'stdRes':modres_data[5]})
    return(MODRES)


def split_MODRES(raw):
    return(  raw.str[6:11].str.replace(' ', ''),  raw.str[12:15].str.replace(' ', ''),     raw.str[16].str.replace(' ', ''), \
            raw.str[17:22].str.replace(' ', ''),     raw.str[22].str.replace(' ', ''),  raw.str[23:27].str.replace(' ', '')
          )


##  Helper function for readjusting residues in MODRES lines
def adjust_MODRES_HETATMs(MODRES, HETATM):
    new_AAs = pd.DataFrame(columns = HETATM.columns)
    for mod_i, mod_row in MODRES.iterrows():
        new_aa = HETATM.copy()
        new_aa = new_aa[new_aa['seqNum'] == mod_row['seqNum']]
        new_aa = new_aa[new_aa['chainID'] == mod_row['chainID']]
        new_aa.loc[:, 'resName'] = mod_row['stdRes']
        new_AAs = new_AAs.append(new_aa)
    return(new_AAs)


def get_MODRES_pdb_df(pdb_file_df):
    ## modifies any HETATM MODRES lines to ATOM lines with stdRes (don't think I need to change SEQRES lines)
    MODRES_DF = get_MODRES_DF(pdb_file_df[pdb_file_df['key'] == 'MODRES'])
    
    ## find index values for atoms in modified residue
    HETATM_DF = get_ATOM_DF(pdb_file_df[pdb_file_df['key'] == 'HETATM'])
    for mod_i, mod_row in MODRES_DF.iterrows():
        ## find index values for atoms in modified residue
        HETATM_MOD_ATOM = HETATM_DF[HETATM_DF['seqNum'] == mod_row['seqNum']]
        HETATM_MOD_ATOM = HETATM_MOD_ATOM[HETATM_MOD_ATOM['chainID'] == mod_row['chainID']]
        HETATM_MOD_ATOM_indexes = HETATM_MOD_ATOM.index.values.tolist()
        ## unmodify ATOM line
        pdb_file_df.loc[HETATM_MOD_ATOM_indexes, 'raw'] = "ATOM  " + pdb_file_df.loc[HETATM_MOD_ATOM_indexes, 'raw'].str[6:17] + str(mod_row['stdRes']) + pdb_file_df.loc[HETATM_MOD_ATOM_indexes, 'raw'].str[20:]
    pdb_file_df['key'] = pdb_file_df['raw'].str[:6].str.replace(' ', '')
    return(pdb_file_df)


## 'ATOM' and 'HETATM' use get_ATOM_DF for a atom level dataframe
def get_SEQADV_DF(SEQADV):
    SEQADV_data = split_SEQADV(SEQADV['raw'])
    SEQADV = pd.DataFrame({
         'idCode':SEQADV_data[0]  , 'resName':SEQADV_data[1],  'chainID':SEQADV_data[2]     \
        ,'seqNum':SEQADV_data[3]  , 'iCode':SEQADV_data[4]  ,  'database':SEQADV_data[5]     \
        ,'dbIdCode':SEQADV_data[6], 'dbRes':SEQADV_data[7]  ,  'dbSeq':SEQADV_data[8]          \
        ,'conflict':SEQADV_data[9]  })
    return(SEQADV)


def split_SEQADV(raw):
    return(  raw.str[6:11].str.replace(' ', ''),  raw.str[12:15].str.replace(' ', ''),     raw.str[16].str.replace(' ', ''), \
            raw.str[17:22].str.replace(' ', ''),     raw.str[23].str.replace(' ', ''),  raw.str[24:27].str.replace(' ', ''), \
            raw.str[29:37].str.replace(' ', ''),  raw.str[39:42].str.replace(' ', ''),  raw.str[43:48].str.replace(' ', ''), \
            raw.str[49:].str.replace(' ', '') )

#########################################################################################################################
#################################################         ATOM_ID       #################################################
#########################################################################################################################
def add_ID_var(ID, var):
    if type(var) == str: ID = ID + '_' + var
    else: ID = ID + '_' + var.astype(str)
    return(ID)

def make_RES_ID(pdb, chainID, resName, seqNum):
    RES_ID = pdb
    RES_ID = add_ID_var(RES_ID, chainID)
    RES_ID = add_ID_var(RES_ID, resName)
    RES_ID = add_ID_var(RES_ID, seqNum)
    return(RES_ID)


def make_ATOM_ID(RES_ID, element, serial):
    ATOM_ID = RES_ID
    ATOM_ID = add_ID_var(ATOM_ID, element)
    ATOM_ID = add_ID_var(ATOM_ID, serial)
    return(ATOM_ID)

def make_tmATOM_ID(RES_ID, serial):
    tmATOM_ID = RES_ID
    tmATOM_ID = add_ID_var(tmATOM_ID, serial)
    return(tmATOM_ID)

def split_ID(ID):
    ## for non series
    if type(ID) == str:
        ID_sp = ID.split('_')
    else:
        ID_sp = ID.str.split('_')
    return(ID_sp)

def pdb_code_from_ATOM_ID(ATOM_ID):
    ATOM_ID_sp = split_ID(ATOM_ID)
    if type(ATOM_ID_sp) == list: pdb_code = ATOM_ID_sp[0]
    else: pdb_code = ATOM_ID_sp.str[0]
    return(pdb_code)

def pdb_code_from_RES_ID(RES_ID):
    RES_ID_sp = split_ID(RES_ID)
    if type(RES_ID_sp) == list: pdb_code = RES_ID_sp[0]
    else: pdb_code = RES_ID_sp.str[0]
    return(pdb_code)

def pdb_name_from_ATOM_ID(ATOM_ID):
    ATOM_ID_sp = split_ID(ATOM_ID)
    if type(ATOM_ID_sp) == list: pdb_name = ATOM_ID_sp[0] + "_" + ATOM_ID_sp[1]
    else: pdb_name = ATOM_ID_sp.str[0] + "_" + ATOM_ID_sp.str[1]
    return(pdb_name)

def chainID_from_ATOM_ID(ATOM_ID):
    ATOM_ID_sp = split_ID(ATOM_ID)
    if type(ATOM_ID_sp) == list: chainID = ATOM_ID_sp[1]
    else: chainID = ATOM_ID_sp.str[1]
    return(chainID)

def chainID_from_RES_ID(RES_ID):
    RES_ID_sp = split_ID(RES_ID)
    if type(RES_ID_sp) == list: chainID = RES_ID_sp[1]
    else: chainID = RES_ID_sp.str[1]
    return(chainID)

def resName_from_ATOM_ID(ATOM_ID):
    ATOM_ID_sp = split_ID(ATOM_ID)
    if type(ATOM_ID_sp) == list: resName = ATOM_ID_sp[2]
    else: resName = ATOM_ID_sp.str[2]
    return(resName)

def resName_from_RES_ID(RES_ID):
    RES_ID_Sp = split_ID(RES_ID)
    if type(RES_ID_Sp) == list: resName = RES_ID_Sp[2]
    else: resName = RES_ID_Sp.str[2]
    return(resName)

def seqNum_from_ATOM_ID(ATOM_ID):
    ATOM_ID_sp = split_ID(ATOM_ID)
    if type(ATOM_ID_sp) == list: seqNum = ATOM_ID_sp[3]
    else: seqNum = ATOM_ID_sp.str[3]
    return(seqNum)

def seqNum_from_RES_ID(RES_ID):
    RES_ID_Sp = split_ID(RES_ID)
    if type(RES_ID_Sp) == list: seqNum = RES_ID_Sp[3]
    else: seqNum = RES_ID_Sp.str[3]
    return(seqNum)

def element_from_ATOM_ID(ATOM_ID):
    ATOM_ID_sp = split_ID(ATOM_ID)
    if type(ATOM_ID_sp) == list: element = ATOM_ID_sp[4]
    else: element = ATOM_ID_sp.str[4]
    return(element)

def serial_from_ATOM_ID(ATOM_ID):
    ATOM_ID_sp = split_ID(ATOM_ID)
    if type(ATOM_ID_sp) == list: serial = ATOM_ID_sp[5]
    else: serial = ATOM_ID_sp.str[5]
    return(serial)

def RES_ID_from_ATOM_ID(ATOM_ID):
    RES_ID = make_RES_ID(pdb_code_from_ATOM_ID(ATOM_ID),
                         chainID_from_ATOM_ID(ATOM_ID),
                         resName_from_ATOM_ID(ATOM_ID),
                         seqNum_from_ATOM_ID(ATOM_ID))
    return(RES_ID)
