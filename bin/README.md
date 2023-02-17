Copyright (C) 2021 University of Kansas

# Contents of MAHOMES-II/bin
This directory contains the non-code files required for saving binary executable models on new machines and evaluating them.
** ML_input/ ** - This folder contains the attempted feature sets during ML model selection.
** MAHOMES_II_sites.csv ** This contains the sites used to train and evaluate MAHOMES II. The columns are as follows
- struc_id: The structure used for that entry. This includes PDB ID and the relax rank. The relax rank is not included for dataset sites since we only train using the top results.
- resName1, resName2, resName3, resName4: The three letter residue code used for metals included in the site	
- seqNum1, seqNum2, seqNum3, seqNum4: The sequence number for metals included in the site
- SITE_ID: An index that can be used to specify the entry. It includes the struc_id information and a # identifier for that site on the given structure. 
- nuclearity: the number of metal atoms included in the site (between 1 and 4)
- Set: ‘data’ for sites in the dataset and ‘test’ for sites in the T-metal-sites10
- pdb_name: The PDB ID and chainID for the entry.
- MAHOMES_match: MAHOMES II redefined sites using the relaxed strucutres (metals within 5 Å). 'True' if the entry contains the same metal atoms as a site used in MAHOMES, which was identified using crystal structures. 'Partial' if the entry is not an exact match to a MAHOMES site, but conatins a metal used in a MAHOMES site. 
- Enzyme: The enzymatic or non-enzymatic target value (updated by Tsite_FPs.xls and kfold_FP.xlsx when relevent)
- meta_site: The handle used for combining predictions from different relaxed structures for the same site.
**  Tsite_FPs.xls and kfold_FP.xlsx ** - Manually examined sites during work on MAHOMES II.
** testset_features1.csv and dataset_features1.csv ** - The feature values used for training and evaluating MAHOMES II. Debugging feature calculations between model selection and publishing MAHOMES II have led to minor differences relative ML_input/.
** AF_sites.csv ** The site identifier's and enzyme/non-enzyme labels for the AlphaFold set.
** AFset_features.csv ** - The feature values for the AlphaFold set sites.


