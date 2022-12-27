Copyright © 2021 University of Kansas

# MAHOMES II
Metal Activity Heuristic of Metalloprotein and Enzymatic Sites (MAHOMES) II - Predicts if a protein bound metal ion is enzymatic or non-enzymatic

## Overview
The ability to distinguish enzyme from non-enzyme sites remains an unsolved, challenging task. We've developed MAHOMES, a machine learning based tool which classifies metals bound to proteins as enzymatic or non-enzymatic. We intend to build on the previous work to make MAHOMES II, a more stable and robust version with a web server.

## System requirements
Feature calculations also require using Rosetta, FindGeo, and bluues which we run using Python 2.7 with CentOS 7 (or red hat now?).

## Installation guide
### set up virtual environment:
```
$ virtualenv --version # check for virtual enviroment
$ pip install virtualenv # download using pip if no version is found
$ virtualenv -p /usr/bin/python3 venv # create new virtual environment
$ source venv/bin/activate # switch to new env 
$ pip install -r requirements.txt # add packages to environment
```
Repeat above proccess for venv2 using python2.7 and requirements2.txt

### Additional required setup
1. Follow FeatureCalculations/README.md to setup third-party feature calculations
2. Update directory paths in setup_ML.sh and driver.sh
3. Train and save MAHOMES II ML models
```
$ ./setup_ML.sh
```
## Instructions for use
To make enzyme and non-enzyme predictions using MAHOMES II
1. Make a new directory and place one or more PDB files in it
2. Run the following command, replacing $JOB_DIR with the directory from step 1
```
. /path-to-repo/MAHOMES-II-server/driver.sh $JOB_DIR
```
The directory should now have predictions.csv as well as all the calculated features.


## References and acknowledgements

MAHOMES II communicates with and/or references the following separate libraries
and packages:
*   [BLUUES](https://doi.org/10.1186/1471-2105-13-S4-S18)
> Fogolari, F., Corazza, A., Yarra, V., Jalaru, A., Viglino, P., & Esposito, G.
     (2012). Bluues: a program for the analysis of the electrostatic
     properties of proteins based on generalized Born radii. BMC
     Bioinformatics, 13(4), S18. doi:10.1186/1471-2105-13-S4-S18
*   [FindGeo](http://metalweb.cerm.unifi.it/tools/findgeo/)
> Andreini, C., Cavallaro, G., & Lorenzini, S. (2012). FindGeo: a tool for
     determining metal coordination geometry. Bioinformatics, 28(12),
     1658-1660. doi:10.1093/bioinformatics/bts246
*   [GHECOM](https://pdbj.org/ghecom/)
> Kawabata, T. (2019). Detection of cave pockets in large molecules: Spaces
     into which internal probes can enter, but external probes from outside
     cannot. Biophysics and physicobiology, 16, 391-406.
     doi:10.2142/biophysico.16.0_391
*   [NumPy](https://numpy.org)
*   [pandas](https://pandas.pydata.org/)
*   [pdb2pqr](https://github.com/Electrostatics/pdb2pqr)
> Jurrus, E., Engel, D., Star, K., Monson, K., Brandi, J., Felberg, L. E., . . . Baker,
     N. A. (2018). Improvements to the APBS biomolecular solvation
     software suite. Protein Science, 27(1), 112-128.
     doi:10.1002/pro.3280
*   [Rosetta software suite](https://www.rosettacommons.org)
> Alford, R. F., Leaver-Fay, A., Jeliazkov, J. R., O’Meara, M. J., DiMaio, F. P., Park,
     H., . . . Gray, J. J. (2017). The Rosetta All-Atom Energy Function for
     Macromolecular Modeling and Design. Journal of Chemical Theory and
     Computation, 13(6), 3031-3048. doi:10.1021/acs.jctc.7b00125
*   [scikit-learn](https://github.com/scikit-learn/scikit-learn)
> Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., Grisel, O., .
     . . Dubourg, V. (2011). Scikit-learn: Machine learning in Python. The
     Journal of Machine Learning Research, 12, 2825-2830.
*   [SciPy](https://scipy.org)


## Contact

The Slusky Lab can be contacted using mahomes@ku.edu for any additional questions about this repo.

## License 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

