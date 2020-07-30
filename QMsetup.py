#!/usr/bin/env python
'''
Program to set up the QM region for CP2K based on an amber topology file.
Requirements: -ambertools (parmed)
              -python 2 or 3 with argparse and pandas library

To run, first make sure you have ambertools in your path:

source /path/to/amber18/amber.sh

Then run
python QMsetup.py -i [topology file] -res [residue index] -link [optional flag to specify linking region] -ligand [name of the ligand in the topology file] -fixLJ [optional flag to fix OH group Lennard Jones charges]

The script will produce 2 output files:
parmed-qmcharge.in  -  parmed script to run to change the charges on the MM region (and optionally fix LJ parameters)
QMresidue-index.inc -  text file that includes index and QM_KIND for atoms in the QM region (and optionally linking atoms) to be included in a CP2K input file

Written by J. McCarty 2020


 This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

from __future__ import print_function
import sys
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
try:
    from parmed.amber import AmberParm
    from parmed.tools import netCharge
    from parmed.tools import printDetails
except ImportError:
    print('make sure you have ambertools in your path:\n\nsource /path/to/amber18/amber.sh\n',file=sys.stderr)
    sys.exit()
import argparse
try:
    import pandas as pd
except ImportError:
    print('python pandas library not found\nInstall with:\n\npip install pandas\n',file=sys.stderr)
    sys.exit()

# function to convert string to list
def Convert(string):
    li = list(string.split(" "))
    return li

# function to get the charges from Amber topology
def get_charges(parm,atomlist):
    size = len(atomlist)
    charges = []
    for i in range(size):
        txt='@'+str(atomlist[i])
        action = netCharge(parm, txt)
        action2str = '%s' % action
        li = Convert(action2str)
        charges.append(float(li[6]))
    return charges

# function to get the difference between two lists
def Diff(list1,list2):
    if(len(list1) > len(list2)):
        return (list(set(list1) - set(list2)))
    else:
        return (list(set(list2) - set(list1)))

# parse arguments from the command line
parser = argparse.ArgumentParser(description='set up charges in QM region')
parser.add_argument('-i',help='amber topology file',type=str, default='complex.prmtop')
parser.add_argument('-res',help='list of residues (example -res 89 77)',type=str,nargs='+',default='')
parser.add_argument('-link',help='add linking residues',action='store_true',default=False)
parser.add_argument('-ligand',help='residue name of the ligand in the topology file (example -ligand CHB)',type=str)
parser.add_argument('-fixLJ',help='optional flag to change LJ parameters for OH groups',action='store_true',default=False)

args = parser.parse_args()

parmfile = args.i
residues = args.res
link = args.link
ligand = args.ligand
write_LJ = args.fixLJ

if(link):
    MM_link = []
    QM_link = []

# print parsed arguments
print('Amber topology file is ', parmfile)
if(len(residues) < 1):
    print('No protein residues will be included in the QM region. Specify residue index to include with -res option.')
else:
    print('Residues to inclued in QM region', residues)
if(link):
    print('Will attempt to get linking atom information')
else:
    print('Will not get linking atom information')

if(ligand):
    print('Will get ligand information for RES',ligand)
else:
    print('No ligand specified. Ligand will be ignored')

if(write_LJ):
    print('Will update LJ parameters for OH groups in the parmed script')
# load parameter file
try:
    parm = AmberParm(parmfile)
except:
    print("Error: amber topology file",parmfile,'could not be opened',file=sys.stderr)
    sys.exit()

# open files for writing
fout_parmed = open('parmed-qmcharge.in','w')
fout_residues = open('QMresidue-index.inc','w')

if(write_LJ):
    fout_parmed.write('changeLJSingleType :WAT@H1 0.3019 0.047\n')
    fout_parmed.write('changeLJSingleType :WAT@H2 0.3019 0.047\n')
    fout_parmed.write('changeLJSingleType :SER@HG 0.3019 0.047\n')
    fout_parmed.write('changeLJSingleType :TYR@HH 0.3019 0.047\n')
    fout_parmed.write('changeLJSingleType :THR@HG1 0.3019 0.047\n')

# an array for all the atom kinds in the QM region
include_atom_kinds = []
total_QM_kinds_index = []

# initialize total QM region charge
total_QM_region_charge = 0

# first get ligand information as needed:
if(ligand):
    txt = ':'+ligand
    residue_info = printDetails(parm,txt)
    residue_info2str = '%s' % residue_info
    residue_data = StringIO(residue_info2str)
    df = pd.read_csv(residue_data, skiprows = 2, sep = ',')
    names=df.columns.tolist()
    name_list=df[names].values.tolist()
    natoms = len(name_list)
    ligand_atoms = []
    atom_kind = []
    atom_name = []
    for j in range(natoms):
        atom_entry = name_list[j][0].split()
        ligand_atoms.append(int(atom_entry[0]))
        name_entry = atom_entry[3]
        atom_name.append(name_entry)
        name_entry = [char for char in name_entry]
        atom_kind.append(name_entry[0])

    if(len(name_list) < 1):
        print("Error: Ligand not found! RES with name ",ligand,'was not found in the topology file',file=sys.stderr)
        sys.exit()

    print('Found ligand', ligand,'! Lignd' , ligand, 'has ', len(name_list),'atom entries with the following indexes to be added to the QM region:')
    print(ligand_atoms)
    for j in range(len(name_list)):
        qm_resnum = ligand_atoms[j]
        ndx = ligand_atoms.index(qm_resnum)
        kind = atom_kind[ndx]
        if(kind not in include_atom_kinds):
            include_atom_kinds.append(kind)
            newrow = []
            newrow.append(qm_resnum)
            total_QM_kinds_index.append(newrow)
        else:
            kind_ndx = include_atom_kinds.index(kind)
            #total_QM_kinds_index[kind_ndx].append(qm_resnum)
            total_QM_kinds_index[kind_ndx].append(qm_resnum)

    # Now deal with charges
    QM_charges = get_charges(parm, ligand_atoms)
    print('total charge on ligand',sum(QM_charges))
    print('An integer charge of',int(round(sum(QM_charges))),'will be added to the QM region')
    total_QM_region_charge += int(round(sum(QM_charges)))
    print('Charge of QM region',total_QM_region_charge)
    print('Done with Ligand')

# loop through residues to include in QM region
for i in range(len(residues)):
    print('Getting information for residue', residues[i],' ',parm.residues[int(residues[i])-1].name)
    txt = ':'+residues[i]
    residue_info = printDetails(parm,txt)
    residue_info2str = '%s' % residue_info
    residue_data = StringIO(residue_info2str)
    df = pd.read_csv(residue_data, skiprows = 2, sep = ',')
    names=df.columns.tolist()
    name_list=df[names].values.tolist()
    natoms = len(name_list)
    total_Atoms = []
    atom_kind = []
    atom_name = []
    for j in range(natoms):
        atom_entry = name_list[j][0].split()
        total_Atoms.append(int(atom_entry[0]))
        name_entry = atom_entry[3]
        atom_name.append(name_entry)
        name_entry = [char for char in name_entry]
        atom_kind.append(name_entry[0])

    print(residue_info2str)
    print('residue', residues[i],' ',parm.residues[int(residues[i])-1].name, 'has ', len(name_list),'atom entries with the following indexes:')
    print(total_Atoms)
    QMdata = raw_input('Select residues in QM region for residue ')
    QMdata = QMdata.replace(','," ")
    QMdata = QMdata.replace('and'," ")
    QM_list=QMdata.split()
    for s in QM_list:
        if(s=="to"):
            ndx=QM_list.index("to")
            lower = int(QM_list[ndx-1])
            upper = int(QM_list[ndx+1])+1
            tmp_list = range(lower,upper)
            QM_list.pop(ndx)
            QM_list.pop(ndx)
            QM_list.pop(ndx-1)
            for j in range(len(tmp_list)):
                QM_list.append(tmp_list[j])

    QM_list = [int(j) for j in QM_list]
    print('Atom index for the QM region')
    print(QM_list)
    check_QM_list = all(elem in total_Atoms for elem in QM_list)

    if check_QM_list:
        print("Pass! all elements of QM list are in the total atom list")
    else:
        print("Error: total atom list does not contain all elements listed in the QM region",file=sys.stderr)
        sys.exit()
    if(len(QM_list)<1):
        print("Error: QM region contains no atoms!",file=sys.stderr)
        sys.exit()
    MM_atoms = Diff(QM_list,total_Atoms)
    print('Atom index for the MM region')
    print(MM_atoms)
    n_MM_atoms = len(MM_atoms)
    n_QM_atoms = len(QM_list)
    for j in range(n_QM_atoms):
        qm_resnum = QM_list[j]
        ndx = total_Atoms.index(qm_resnum)
        kind = atom_kind[ndx]
        if(kind not in include_atom_kinds):
            include_atom_kinds.append(kind)
            newrow = []
            newrow.append(qm_resnum)
            total_QM_kinds_index.append(newrow)
        else:
            kind_ndx = include_atom_kinds.index(kind)
            #total_QM_kinds_index[kind_ndx].append(qm_resnum)
            total_QM_kinds_index[kind_ndx].append(qm_resnum)

    if(link):
        found_MM_CA = False
        found_QM_CB = False
        if('CA' in atom_name):
            trial_MM_ndx = atom_name.index('CA')
            trial_MM = total_Atoms[trial_MM_ndx]
            if(trial_MM in MM_atoms):
                found_MM_CA = True
        if(found_MM_CA):
            txt = 'Select index number for the MM linking atom in residue '+str(residues[i])+' '+str(parm.residues[int(residues[i])-1].name)+'   (Recommended CA atom index'+' '+str(trial_MM)+')  '
        else:
            print('Note: no CA atom was found in the MM region')
            txt = 'Select index number for the MM linking atom in residue '+str(residues[i])+' '+str(parm.residues[int(residues[i])-1].name)+'   '
        MMlinkingatom = raw_input(txt)
        if('CB' in atom_name):
            trial_QM_ndx = atom_name.index('CB')
            trial_QM = total_Atoms[trial_QM_ndx]
            if(trial_QM in QM_list):
                found_QM_CB = True
        if(found_QM_CB):
            txt = 'Select index number for the QM linking atom in residue '+str(residues[i])+' '+str(parm.residues[int(residues[i])-1].name)+'   (Recommended CB atom index'+' '+str(trial_QM)+')  '
        else:
            print('Note: no CB atom was found in the QM region')
            txt = 'Select index number for the QM linking atom in residue '+str(residues[i])+' '+str(parm.residues[int(residues[i])-1].name)+'   '
        QMlinkingatom = raw_input(txt)
        if(len(QMlinkingatom)<1 or len(MMlinkingatom)<1):
            print('One or more linking atoms was not specified. Will not write linking atoms for this residue.')
        else:
            MM_index=total_Atoms.index(int(MMlinkingatom))
            QM_index=total_Atoms.index(int(QMlinkingatom))
            print('Will link atoms ',MMlinkingatom,' ',atom_name[MM_index],' with ',QMlinkingatom,' ',atom_name[QM_index])
            MM_link.append(MMlinkingatom)
            QM_link.append(QMlinkingatom)
            if(int(QMlinkingatom) not in QM_list):
                print("Error: ",QMlinkingatom, ' is not listed in the QM list',file=sys.stderr)
                sys.exit()
            if(int(MMlinkingatom) not in MM_atoms):
                print("Error: ",MMlinkingatom, ' is not listed in the MM list',file=sys.stderr)
                sys.exit()

    # Now deal with charges
    QM_charges = get_charges(parm, QM_list)
    total_QM_region_charge += int(round(sum(QM_charges)))
    total_charges = get_charges(parm, total_Atoms)
    MM_charges = get_charges(parm, MM_atoms)
    residual_charge = sum(MM_charges)
    print('Now gathering charge information')
    print('sum of charges on QM groups', sum(QM_charges))
    print('An integer charge of',int(round(sum(QM_charges))),'will be added to the QM region')
    print('sum of residual charge on MM groups', residual_charge)
    print('Sanity Check: total QM group charges + total group MM charges ',sum(QM_charges)+residual_charge," and total charge on residue",residues[i],' ',parm.residues[int(residues[i])-1].name,' ', sum(total_charges))
    if(sum(QM_charges)+residual_charge - sum(total_charges) < 1e-8):
        print('Pass sanity check')
    else:
        print('Warning: Sum of QM group charges + MM group charges does not equal the total residue charge')
    print('number of atoms in MM region',n_MM_atoms)
    print('Will redistribute residual charge of ',residual_charge,' over ', n_MM_atoms, ' atoms in the MM region')

    add_charge = -1.0*float(residual_charge)/(float(n_MM_atoms))
    print('Adding charge ', add_charge, 'to each MM atom')

    new_MM_charge = 0.0
    for j in range(n_MM_atoms):
        txt='@'+str(MM_atoms[j])
        current_charge = get_charges(parm,[MM_atoms[j]])
        print('current charge on MM residue',MM_atoms[j],' ',current_charge[0])
        new_charge = current_charge[0] + add_charge
        new_MM_charge += new_charge
        wtxt = 'change charge '+txt+' '+str(new_charge)+'\n'
        fout_parmed.write(wtxt)
        print(wtxt)

    print('new residual MM charge will be', new_MM_charge)
    print('Done with ',residues[i],' ',parm.residues[int(residues[i])-1].name)
print('Done with all residues')

print('Atom types to be included in the QM region',include_atom_kinds)
print('Total integer charge on the QM region', total_QM_region_charge)
print('################################################################')
print('                     1.   RUN PARMED          ')
print('################################################################')
print('Run parmed script with:    ')
print('\n')
print('          parmed ',parmfile,'-i parmed-qmcharge.in')
print('\n')
print('################################################################')
print('    2.  MAKE THE FOLLOWING CHANGES TO YOUR CP2K INPUT FILE        ')
print('################################################################')
print('Set CHARGE ', total_QM_region_charge, 'in &FORCE_EVAL / &DFT section of CP2K')
print('Set the PARM_FILE_NAME system_qm-charge.parm7 in your CP2K input file')
print('Set the CONN_FILE_NAME system_qm-charge.parm7 in your CP2K input file')
print('Copy the contents of QMresidue-index.inc to the &FORCE_EVAL / &QMMM / QM_KIND section of your CP2K input file')
print('or include the following in under the &QMMM section:')
print("@INCLUDE './QMresidue-index.inc'")
for k in range(len(include_atom_kinds)):
    txt = '&QM_KIND  '+include_atom_kinds[k]+'\n'
    fout_residues.write(txt)
    for j in range(len(total_QM_kinds_index[k])):
        txt = 'MM_INDEX '+str(total_QM_kinds_index[k][j])+'\n'
        fout_residues.write(txt)
    fout_residues.write('&END QM_KIND\n')

if(link):
    for i in range(len(MM_link)):
        fout_residues.write('&LINK\n')
        txt = 'MM_INDEX  '+str(MM_link[i])+'\n'
        fout_residues.write(txt)
        txt = 'QM_INDEX  '+str(QM_link[i])+'\n'
        fout_residues.write(txt)
        fout_residues.write('LINK_TYPE IMOMM\n')
        fout_residues.write('&END LINK\n')

fout_parmed.write('outparm system_qm-charge.parm7\n')
fout_parmed.write('quit\n')

fout_parmed.close()
fout_residues.close()
