#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Extract Dynamical matrix from ABINIT output file
#
# Author; Eugene Roginskii
#

atom_data = [ 
    [  0, "X", "X", 0], # 0
    [  1, "H", "Hydrogen", 1.00794], # 1
    [  2, "He", "Helium", 4.002602], # 2
    [  3, "Li", "Lithium", 6.941], # 3
    [  4, "Be", "Beryllium", 9.012182], # 4
    [  5, "B", "Boron", 10.811], # 5
    [  6, "C", "Carbon", 12.0107], # 6
    [  7, "N", "Nitrogen", 14.0067], # 7
    [  8, "O", "Oxygen", 15.9994], # 8
    [  9, "F", "Fluorine", 18.9984032], # 9
    [ 10, "Ne", "Neon", 20.1797], # 10
    [ 11, "Na", "Sodium", 22.98976928], # 11
    [ 12, "Mg", "Magnesium", 24.3050], # 12
    [ 13, "Al", "Aluminium", 26.9815386], # 13
    [ 14, "Si", "Silicon", 28.0855], # 14
    [ 15, "P", "Phosphorus", 30.973762], # 15
    [ 16, "S", "Sulfur", 32.065], # 16
    [ 17, "Cl", "Chlorine", 35.453], # 17
    [ 18, "Ar", "Argon", 39.948], # 18
    [ 19, "K", "Potassium", 39.0983], # 19
    [ 20, "Ca", "Calcium", 40.078], # 20
    [ 21, "Sc", "Scandium", 44.955912], # 21
    [ 22, "Ti", "Titanium", 47.867], # 22
    [ 23, "V", "Vanadium", 50.9415], # 23
    [ 24, "Cr", "Chromium", 51.9961], # 24
    [ 25, "Mn", "Manganese", 54.938045], # 25
    [ 26, "Fe", "Iron", 55.845], # 26
    [ 27, "Co", "Cobalt", 58.933195], # 27
    [ 28, "Ni", "Nickel", 58.6934], # 28
    [ 29, "Cu", "Copper", 63.546], # 29
    [ 30, "Zn", "Zinc", 65.38], # 30
    [ 31, "Ga", "Gallium", 69.723], # 31
    [ 32, "Ge", "Germanium", 72.64], # 32
    [ 33, "As", "Arsenic", 74.92160], # 33
    [ 34, "Se", "Selenium", 78.96], # 34
    [ 35, "Br", "Bromine", 79.904], # 35
    [ 36, "Kr", "Krypton", 83.798], # 36
    [ 37, "Rb", "Rubidium", 85.4678], # 37
    [ 38, "Sr", "Strontium", 87.62], # 38
    [ 39, "Y", "Yttrium", 88.90585], # 39
    [ 40, "Zr", "Zirconium", 91.224], # 40
    [ 41, "Nb", "Niobium", 92.90638], # 41
    [ 42, "Mo", "Molybdenum", 95.96], # 42
    [ 43, "Tc", "Technetium", 0], # 43
    [ 44, "Ru", "Ruthenium", 101.07], # 44
    [ 45, "Rh", "Rhodium", 102.90550], # 45
    [ 46, "Pd", "Palladium", 106.42], # 46
    [ 47, "Ag", "Silver", 107.8682], # 47
    [ 48, "Cd", "Cadmium", 112.411], # 48
    [ 49, "In", "Indium", 114.818], # 49
    [ 50, "Sn", "Tin", 118.710], # 50
    [ 51, "Sb", "Antimony", 121.760], # 51
    [ 52, "Te", "Tellurium", 127.60], # 52
    [ 53, "I", "Iodine", 126.90447], # 53
    [ 54, "Xe", "Xenon", 131.293], # 54
    [ 55, "Cs", "Caesium", 132.9054519], # 55
    [ 56, "Ba", "Barium", 137.327], # 56
    [ 57, "La", "Lanthanum", 138.90547], # 57
    [ 58, "Ce", "Cerium", 140.116], # 58
    [ 59, "Pr", "Praseodymium", 140.90765], # 59
    [ 60, "Nd", "Neodymium", 144.242], # 60
    [ 61, "Pm", "Promethium", 0], # 61
    [ 62, "Sm", "Samarium", 150.36], # 62
    [ 63, "Eu", "Europium", 151.964], # 63
    [ 64, "Gd", "Gadolinium", 157.25], # 64
    [ 65, "Tb", "Terbium", 158.92535], # 65
    [ 66, "Dy", "Dysprosium", 162.500], # 66
    [ 67, "Ho", "Holmium", 164.93032], # 67
    [ 68, "Er", "Erbium", 167.259], # 68
    [ 69, "Tm", "Thulium", 168.93421], # 69
    [ 70, "Yb", "Ytterbium", 173.054], # 70
    [ 71, "Lu", "Lutetium", 174.9668], # 71
    [ 72, "Hf", "Hafnium", 178.49], # 72
    [ 73, "Ta", "Tantalum", 180.94788], # 73
    [ 74, "W", "Tungsten", 183.84], # 74
    [ 75, "Re", "Rhenium", 186.207], # 75
    [ 76, "Os", "Osmium", 190.23], # 76
    [ 77, "Ir", "Iridium", 192.217], # 77
    [ 78, "Pt", "Platinum", 195.084], # 78
    [ 79, "Au", "Gold", 196.966569], # 79
    [ 80, "Hg", "Mercury", 200.59], # 80
    [ 81, "Tl", "Thallium", 204.3833], # 81
    [ 82, "Pb", "Lead", 207.2], # 82
    [ 83, "Bi", "Bismuth", 208.98040], # 83
    [ 84, "Po", "Polonium", 0], # 84
    [ 85, "At", "Astatine", 0], # 85
    [ 86, "Rn", "Radon", 0], # 86
    [ 87, "Fr", "Francium", 0], # 87
    [ 88, "Ra", "Radium", 0], # 88
    [ 89, "Ac", "Actinium", 0], # 89
    [ 90, "Th", "Thorium", 232.03806], # 90
    [ 91, "Pa", "Protactinium", 231.03588], # 91
    [ 92, "U", "Uranium", 238.02891], # 92
    [ 93, "Np", "Neptunium", 0], # 93
    [ 94, "Pu", "Plutonium", 0], # 94
    [ 95, "Am", "Americium", 0], # 95
    [ 96, "Cm", "Curium", 0], # 96
    [ 97, "Bk", "Berkelium", 0], # 97
    [ 98, "Cf", "Californium", 0], # 98
    [ 99, "Es", "Einsteinium", 0], # 99
    [100, "Fm", "Fermium", 0], # 100
    [101, "Md", "Mendelevium", 0], # 101
    [102, "No", "Nobelium", 0], # 102
    [103, "Lr", "Lawrencium", 0], # 103
    [104, "Rf", "Rutherfordium", 0], # 104
    [105, "Db", "Dubnium", 0], # 105
    [106, "Sg", "Seaborgium", 0], # 106
    [107, "Bh", "Bohrium", 0], # 107
    [108, "Hs", "Hassium", 0], # 108
    [109, "Mt", "Meitnerium", 0], # 109
    [110, "Ds", "Darmstadtium", 0], # 110
    [111, "Rg", "Roentgenium", 0], # 111
    [112, "Cn", "Copernicium", 0], # 112
    [113, "Uut", "Ununtrium", 0], # 113
    [114, "Uuq", "Ununquadium", 0], # 114
    [115, "Uup", "Ununpentium", 0], # 115
    [116, "Uuh", "Ununhexium", 0], # 116
    [117, "Uus", "Ununseptium", 0], # 117
    [118, "Uuo", "Ununoctium", 0], # 118
    ]


def chunks(l, n):
    if n < 1:
        n = 1
    return [l[i:i + n] for i in range(0, len(l), n)]

def strType(var):
    try:
        if int(var) == float(var):
            return 'int'
    except:
        try:
            float(var)
            return 'float'
        except:
            return 'str'

def isNumeric(var):
    isnum=strType(var)
    if (isnum == 'int' or isnum == 'float'):
        return (True)
    return(False)



import numpy as np
from math import sqrt
import sys
import re
import argparse

# To produce Dynamic matrix with capabillity of PHONOPY code one need following constants
Ha2Ev=27.2107
Angstr2Bohr=1.8897261245650618

parser = argparse.ArgumentParser(description='Script to extract Dynamical matrix from abinit output file')

parser.add_argument("-i", "--input", action="store", type=str, dest="abinit_fn",  help="Abinit output filename")
parser.add_argument("-o", "--output", action="store", type=str, dest="out_fn",  default='qpoints.yaml', help="Output filename in yaml format")
parser.add_argument("-d", "--dset", action="store", type=int, dest="dset", default=0, help="Dataset number to extract dynmat from")

args = parser.parse_args()

if (args.abinit_fn == None):
    print('Error. No Abinit input filename was given.')
    sys.exit(1)

try:
    abinit_fh = open(args.abinit_fn, 'r')
except IOError:
    print "ERROR Couldn't open abinit file, exiting...\n"
    sys.exit(1)

for line in abinit_fh:
    if 'Echo of variables that govern the present computation' in line:
        break
natom=0
for line in abinit_fh:
    if 'natom' in line:
        if(isNumeric(line.split()[1])):
            natom=int(line.split()[1])
            break
        else:
            print('Error getting number of atoms value')
            sys.exit(1)

if (natom==0):
    print('Error getting number of atoms value')
    sys.exit(1)

abinit_fh.seek(0)
for line in abinit_fh:
    if 'Echo of variables that govern the present computation' in line:
        break

#Python bug fix. Just call seek once without offset
abinit_fh.seek(0,1)

while True:
    line=abinit_fh.readline()
    if not line: break
    if (re.match('^\s*#',line)):
        continue
    if '#' in line:
        line=re.sub('\s*#.*','',line)

    if ('dataset' in line.lower()):
        if (re.match('^=+\s*DATASET\s+(\d+).*',line)):
            break

    if 'ntypat' in line:
         ntypat=int(line.split()[1])
    elif 'typat' in line:
        if(re.match('^\s*typat\s+\d',line) is None):
            line=abinit_fh.readline()
            if(re.match('^\s*\d',line) is None): 
                print('Error reading typat variable in abinit input file')
                sys.exit(1)
            else:
                typat=[int(line.split()[i]) for i in range(natom)]
        else:
            typat=[int(line.split()[i+1]) for i in range(natom)] 
    elif 'znucl' in line:
        znuclstr=line.split()[1:]
        print(znuclstr)
#Yep abinit print charge as a float
        znucl=[int(float(znuclstr[i])) for i in range(len(znuclstr)) ]
        print('znucl= %s' % znucl)


print('typat: %s' % typat)

# This is how to get mass: here N is the type of atom
#print(atom_data[znucl[typat[N]-1]][3])

#sys.exit(1)

dm=[]
abinit_fh.seek(0)
if(args.dset!=0):
#Skip header
    for line in abinit_fh:
        if 'Echo of variables that govern the present computation' in line:
            break

    for line in abinit_fh:
        if ('dataset' in line.lower()):
            if (re.match('^=+\s*DATASET\s+(\d+).*',line)):
                if (isNumeric(re.match('^=+\s*DATASET\s+(\d+).*',line).group(1))):
                    if(int(re.match('^=+\s*DATASET\s+(\d+).*',line).group(1)) == args.dset):
                        print('Extracting dynmat from DATASET %d' % args.dset)
                        break


for line in abinit_fh:
    if ('dataset' in line.lower()):
            if (re.match('^=+\s*DATASET\s+(\d+).*',line)):
                print('No dynmat in DATASET %d' % args.dset)
                sys.exit(1)
    if 'dynamical matrix, in cartesian coordinates' in line.lower():
        break

row=0
col=0
dmline=[]
mu=0
for line in abinit_fh:
    if (len(line.split())<5):
        continue
    elif ((isNumeric(line.split()[4])!=True) or (isNumeric(line.split()[5]) != True) or
                    (isNumeric(line.split()[1])!=True) or (isNumeric(line.split()[3]) != True)):
        continue
# atom_data[znucl[typat[N]-1]][3]
    mu=atom_data[znucl[typat[int(line.split()[1])-1]-1]][3] * atom_data[znucl[typat[int(line.split()[3])-1]-1]][3]
    if (mu == 0):
        print('Error. Mass is zero')
        sys.exit(1)
    dmline.append([(float(line.split()[4])/sqrt(mu)*Ha2Ev*Angstr2Bohr), (float(line.split()[5])/sqrt(mu)*Ha2Ev*Angstr2Bohr)])
    row=row+1
    if (row>=natom*3):
        row=0
        col=col+1
        dm.append(dmline)
        dmline=[]
    if(col>=natom*3):
        break

if (len(dm)<3*natom):
    print('Error reading dynmat:')

try:
    out_fh = open(args.out_fn, 'w')
except IOError:
    print "ERROR Couldn't open output file for writing, exiting...\n"
    sys.exit(1)

out_fh.write('natom: %d\n' % natom)
out_fh.write('phonon:\n    dynamical_matrix:\n')

for i in range(natom*3):
    out_fh.write('    -  [')
    out_fh.write(', '.join('% 12.10f, % 12.10f' % (d[0],d[1]) for d in dm[i]))
    out_fh.write(']\n')
