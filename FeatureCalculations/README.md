# MAHOMES II physiochemical feature calculations


## Overview


## System requirements
Feature calculations also require using Rosetta, FindGeo, and bluues which we run using Python 2.7 with CentOS 7 (or red hat now?).

## set-up guide
### Successfuly install the following third party software tools:
1. Rosetta 3.13
    - Download:  https://www.rosettacommons.org/software/license-and-download
    - Build instructions: https://new.rosettacommons.org/docs/latest/build_documentation/Build-Documentation 
2. FindGeo: python version can be found at http://metalweb.cerm.unifi.it/tools/findgeo/
3. bluues: can be downloaded in the supplementary from its publication ( https://doi.org/10.1186/1471-2105-13-S4-S18 )
4. pdb2pqr (should already be set-up from venv)

note: FindGeo and bluues binaries may require additional dependencies depending on the system being used.

### change something
1. add the following to your bashrc or bash_profile
export ROSETTA3="/path/to/rosetta/main/source"
export BLUUES_DIR="/path/to/bluues"
export GHECOM_DIR="/path/to/ghecom/"
export MAHOMES_II_DIR="/path/to/MAHOMES_II"

2. change FINDGEO_DIR in ThirdPartyCalculations/findgeo.py to point to the correct location on your machine.

3. run ThirdPartyCalculations/testRun.sh check if ThirdPartyCalculations/check_setup/output/ is correct.

 
