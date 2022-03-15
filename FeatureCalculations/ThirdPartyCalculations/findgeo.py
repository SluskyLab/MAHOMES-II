#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2010-2011
# findgeo is distributed under the terms of the GNU General Public License
#
"""
    findgeo is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    findgeo is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import os
import shutil
import commands
import subprocess
from optparse import OptionParser
import urllib
import numpy as np

#FINDGEO_DIR="/panfs/pfs.local/work/slusky/CommonPrograms/findgeo"
FINDGEO_DIR="/var/www/apps/mahomes/findgeo"

## p3d should be located in same directory as findgeo binary
sys.path.insert(0, "%s" % FINDGEO_DIR)
from p3d.protein import Protein

# This script analyzes the PDB coordinate files and extracts sites
def main(code, input_file, output_dir, distance, metal, pdbID):

    ref_metalList = ['LI', 'BE', 'NA', 'MG', 'AL', 'K', 'CA', 'SC', 'TI', 'V', 'CR', 'MN',
            'FE', 'CO', 'NI', 'CU', 'ZN', 'GA', 'GE', 'RB', 'SR', 'Y', 'ZR', 'NB',
            'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'CS', 'BA',
            'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER',
            'TM', 'YB', 'LU', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG',
            'TL', 'PB', 'BI', 'PO', 'FR', 'RA', 'AC', 'TH', 'PA', 'U', 'NP', 'PU',
            'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO', 'LR', 'RF', 'DB', 'SG',
            'AS']

    """Main function. Prepare the find geometry steps."""
    if input_file:
        code = os.path.join('%s' % code)
        if not os.path.isfile(code):
            sys.exit('%s: file not found' % code)

    #read in pdb file and convert to p3d object
    try:        
        pdb = Protein(code)
    except:
        sys.exit('PDB in bad format.')

    if metal:
        metalList = []
        metalList.append(metal)
    else:
        metalList = ref_metalList

    metals = find_metals(pdb, metalList)

    if metals:
        out = open(os.path.join(output_dir,'%s.findgeo'%pdbID), 'w')
        ligands = find_ligands(metals, pdb, metalList, distance, not_donors)
        for key in ligands:
            print "%s_%s%s_%s_%s --> %d ligands found." %(key.atype.strip(), key.resid, key.altConf, key.idx, key.chain, len(ligands[key]))

            if len(ligands[key]) >= 2 and len(ligands[key]) <= 9:
                #print "Determining coordination geometry..."
                pdb_name = os.path.join(output_dir,'%s_%s_%s_%s.pdb' % (key.atype.strip(), key.resid, key.idx, key.chain))
                outfile = os.path.join(output_dir, pdb_name)
                temp = open(outfile, 'w')
                for line in pdb.query(key):
                    temp.writelines('%s\n' % line.output())
                for element in ligands[key]:
                    temp.writelines('%s\n' % element.output())
                temp.close()
                find_geo_input(pdb_name, output_dir, out, pdbID, ligands)
            elif len(ligands[key]) < 2:
                print "Too few donor atoms identified. At least 2 donor atoms must be present."
            else:
                print "Too many donor atoms identified. The maximum allowed coordination number is 9."
        out.close()
    else:
        if metal:
            sys.exit('%s not found in the PDB file' % metal)
        sys.exit('No metals found in the PDB file')

def find_metals(pdb, metalList):
    """Find metals in the pdb"""
    metals = [atom for atom in pdb.query('non-protein and model 1') if atom.elementType.upper().strip() in metalList]
    return metals

def find_ligands(metals, pdb, metalList, maxDist, not_donors):
    """Find ligands in the pdb"""
    ligands = {}

    for metal in metals:
        ligands[metal] = []
        temp = pdb.query('model 1 & within {0} of'.format(maxDist), metal)

        for atom in temp:
            if atom.elementType.upper().strip() not in not_donors and atom.elementType.upper().strip() not in metalList:
                ligands[metal].append(atom)

        if not ligands[metal]:
            print "No ligands found for metal %s_%s_%s_%s (try to change distance threshold (-t) or excluded donors (-e))" % (metal.atype.strip(), metal.resid, metal.idx, metal.chain)
    return ligands

def find_geo_input(pdb_name, output_dir, out, pdbID, these_ligands):
    """Find geometries and print the best"""

    geometries = {'lin': 'linear',
        'irr': 'irregular',
        'trv': 'trigonal plane with a vacancy',
        'tri': 'trigonal plane',
        'tev': 'tetrahedron with a vacancy',
        'spv': 'square plane with a vacancy',
        'tet': 'tetrahedron',
        'spl': 'square plane',
        'bva': 'trigonal bipyramid with a vacancy (axial)',
        'bvp': 'trigonal bipyramid with a vacancy (equatorial)',
        'pyv': 'square pyramid with a vacancy (equatorial)',
        'spy': 'square pyramid',
        'tbp': 'trigonal bipyramid',
        'tpv': 'trigonal prism with a vacancy',
        'oct': 'octahedron',
        'tpr': 'trigonal prism',
        'pva': 'pentagonal bipyramid with a vacancy (axial)',
        'pvp': 'pentagonal bipyramid with a vacancy (equatorial)',
        'cof': 'octahedron, face monocapped with a vacancy (capped face)',
        'con': 'octahedron, face monocapped with a vacancy (non-capped face)',
        'ctf': 'trigonal prism, square-face monocapped with a vacancy (capped face)',
        'ctn': 'trigonal prism, square-face monocapped with a vacancy (non-capped edge)',
        'pbp': 'pentagonal bipyramid',
        'coc': 'octahedron, face monocapped',
        'ctp': 'trigonal prism, square-face monocapped',
        'hva': 'hexagonal bipyramid with a vacancy (axial)',
        'hvp': 'hexagonal bipyramid with a vacancy (equatorial)',
        'cuv': 'cube with a vacancy',
        'sav': 'square antiprism with a vacancy',
        'hbp': 'hexagonal bipyramid',
        'cub': 'cube',
        'sqa': 'square antiprism',
        'boc': 'octahedron, trans-bicapped',
        'bts': 'trigonal prism, square-face bicapped',
        'btt': 'trigonal prism, triangular-face bicapped',
        'ttp': 'trigonal prism, square-face tricapped',
        'csa': 'square antiprism, square-face monocapped'}

    dir_name = pdb_name.split('.pdb')[0]
    temp = dir_name.split('/')
    name = temp[len(temp)-1]
    if os.path.exists('%s' % dir_name):
        if overwrite:
            try:
                shutil.rmtree('%s' % dir_name)
            except:
                sys.exit('Cannot remove the existing directory %s: check permissions' % dir_name)
        else:
            sys.exit('%s: is an existing directory... Remove/rename it or use the -o option' % dir_name)

    try:
        os.mkdir(dir_name)
    except:
        sys.exit('Cannot create the directory %s: check permissions' % dir_name)

    shutil.move(pdb_name, os.path.join(dir_name, 'findgeo.input'))
    cmd = '%s/findgeo %s' % (FINDGEO_DIR, dir_name)
    geom = commands.getoutput(cmd).split()

    ## CUSTOM OUTPUT ##
    geom_out = geom[0]
    rmsd = 100
    regular_out = "Irr"

    if geom_out != 'irr':
        rmsd = subprocess.check_output(["grep", "^%s"%geom_out, "%s/findgeo.out"%dir_name])
        rmsd = rmsd.strip().split("\n")[0].split()[-1]
        regular_out = geom[1][1:-1]
        #print 'Best fit geometry: %s %s RMSD:%s' % (geometries[geom[0]], geom[1], rmsd)
    #else:
        #print 'Irregular geometry'
    print "attempting to find %s"%name
    ligand_cnt = np.zeros(4) #N, O, S, other
    for key in these_ligands:
       cur_name="%s_%s_%s_%s" %(key.atype.strip(), key.resid, key.idx, key.chain)
       print "\t%s"%cur_name

       if cur_name == name:
           for element in these_ligands[key]:
               atom_id = element.atype
               if "N" in atom_id:
                   ligand_cnt[0] += 1
               elif "O" in atom_id:
                   ligand_cnt[1] += 1
               elif "S" in atom_id:
                   ligand_cnt[2] += 1
               else:
                   ligand_cnt[3] += 1

    ligand_out="\t".join(map(str, ligand_cnt.tolist()))
    out.write('%s\t%s\t%s\t%s\t%s\n' % (name, geom_out, regular_out, rmsd, ligand_out))


if __name__ == '__main__':

    usage = "usage: %prog -p pdbfile [-c pdbcode] [-t threshold] [-m metal] [-e excluded_donors] [-w workdir] [-o overwrite] [-i pdb_id]"

    parser = OptionParser(usage)
    parser.add_option("-w", "--wdir", dest="wd",
            help="Directory where to find or download the input PDB file and to write outputs. Default is ./",
            metavar="WORKDIR")
    parser.add_option("-p", "--pdb_file", dest="pdb_file",
            help="Local input PDB file.",
            metavar="FILE")
    parser.add_option("-c", "--pdb_code", dest="pdb_code",
            help="PDB code of input PDB file to be downloaded from the web.",
            metavar="PDBCODE")
    parser.add_option("-t", "--threshold", dest="threshold",
            help="Coordination distance threshold. Default is 2.8 A.",
            metavar="THRESHOLD")
    parser.add_option("-m", "--metal", dest="metal",
            help="Chemical symbol of the metal of interest. Default is all metals.",
            metavar="METAL")
    parser.add_option("-e", "--excluded_donors", dest="notdonors",
            help="Chemical symbols of the atoms (separated by commas) excluded from metal ligands. Default is C and H.",
            metavar="NOTDONORS")
    parser.add_option("-o", "--overwrite", dest="overwrite",
            action="store_true", default=False,
            help="Overwrite existing files and directories.",
            metavar="OVERWRITE")
    parser.add_option("-i", "--pdb_id", dest="pdb_id",
            help="PDB code of files to save summaries to.",
            metavar="PDBID")

    (options, args) = parser.parse_args()

    if not options.pdb_file and not options.pdb_code:
        parser.error("A PDB file or PDB code is required.")

    if  options.pdb_file and options.pdb_code:
        parser.error("A PDB file or PDB code is required, not both.")

    if options.wd:
        wd = os.path.abspath(options.wd) + '/'
    else:
        wd = './'

    if options.metal:
        metal = options.metal.upper()
        if metal not in ref_metalList:
            parser.error('%s: invalid metal.' % metal)
    else:
        metal = None

    if options.threshold:
        try:
            threshold = float(options.threshold)
        except ValueError:
            sys.exit('Invalid threshold. This must be a number.')
    else:
        threshold = 2.8

    if options.pdb_file:
        pdb = options.pdb_file
        input_file = True
    elif options.pdb_code:
        pdb = options.pdb_code
        input_file = False
        if len(pdb) != 4:
            parser.error('Invalid PDB code. Must be four characters.')

    if options.pdb_id:
        pdb_id=options.pdb_id
    else:
        pdb_id="FindGeo"

    if options.notdonors:
        not_donors = options.notdonors.replace(' ','').split(',')
    else:
        not_donors = ['C', 'H']

    overwrite = options.overwrite

    main(pdb, input_file, wd, threshold, metal, pdb_id)
    
