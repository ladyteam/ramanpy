#!/usr/bin/python3
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
# The raman tensor calculation based on finite displacement method
# Please cite the following paper:
# DOI: 10.1021/acs.jpcc.5b05540
#
# The code is based on PHONOPY API https://phonopy.github.io/phonopy/
#
# Author: Eugene Roginskii

import phonopy
from phonopy.interface.castep import write_castep
from phonopy.interface.vasp import write_vasp
from phonopy.interface.abinit import write_abinit
from phonopy.structure.atoms import PhonopyAtoms as Atoms
import argparse
import numpy as np
import sys
from math import sqrt
import re
from math import (
        sqrt, pi, exp)
import os
import xml.etree.cElementTree as etree


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

def chunks(l, n):
    if n < 1:
        n = 1
    return [l[i:i + n] for i in range(0, len(l), n)]


def direct2cart(directpos,basis):
    cart=[]
    for atdirect in directpos:
        cart.append(np.dot(np.transpose(basis),atdirect))
    return np.array(cart)

def cart2direct(cart,basis):
    direct=[]
    for atcart in cart:
        direct.append(np.dot(np.linalg.inv(np.transpose(basis)),atcart))
    return np.array(direct)

def _iterparse(fname, tag=None):
    for event, elem in etree.iterparse(fname):
            if tag is None or elem.tag == tag:
                        yield event, elem

def gencastep(fn,species,basis,cart, magmoms=None, commentstr='', ir=''):
    cell = Atoms(cell=basis, symbols=species, 
            scaled_positions=cart2direct(cart,basis),magmoms=magmoms)
    try:
        write_castep(fn, cell)
    except:
        print('Error writing file %s' % fn)
    try:
        with open(fn,"r") as fh:
            lines = fh.readlines() #read
        with open(fn, "w") as fh:
            fh.write('#%s; IR: %s\n' % (commentstr,ir))
            fh.writelines(lines) #write back
    except IOError:
        print("ERROR adding comment to output file %s" % fn)
        sys.exit(1)

####### Extract permittivity from DFPT #######
def get_epsilon_dfpt_castep(fn):
# Optical Permittivity (f->infinity)
    try:
        fh = open(fn, 'r')
    except IOError:
        print("ERROR. Couldn't open CASTEP output file")
        return(-1)
    epsilon=[0 for i in range (9)]
    for line in fh:
        if("Optical Permittivity" in line and  "f->infinity" in line):
            break
    fh.readline()
    for i in range(3):
        line=fh.readline()
        for j in range(3):
            epsilon[i*3+j]=float(line.split()[j])
    print(epsilon)
    return(np.array(epsilon))

####### Extract permittivity from OPTIC #######
def get_epsilon_optic_castep(fn):
# # Component            1
    try:
        fh = open(fn, 'r')
    except IOError:
        print("ERROR. Couldn't open CASTEP output file")
        return(-1)
    e=[0 for i in range (6)]
    idx=1
    for i in range(6):
        for line in fh:
            if ('Component' in line):
                idx=int(line.lower().split('component')[1])-1
                break
        for line in fh:
            if (len(line.split())>=3):
                e[idx]=float(line.split()[1])
                break

#      6-components:   0  1  2  3  4  5 
#                     xx yy zz xy xz yz
#      9-components:  xx   xy   xz   yx   yy   yz   zx   zy   zz
    return(np.array([e[0],e[3],e[4],e[3],e[1],e[5],e[4],e[5],e[2]]))

def genposcar(fn,modenum,freq,basis,species,cart,ir=''):
    cell = Atoms(cell=basis, symbols=species, 
            scaled_positions=cart2direct(cart,basis))
    try:
        write_vasp(fn, cell)
    except:
        print('Error writing file %s' % fn)
    try:
        with open(fn,"r") as fh:
            lines = fh.readlines() #read
        with open(fn, "w") as fh:
            fh.write('%s\n' % 'Epsilon calculation for %d mode, %f cm-1 IR: %s' % (modenum,freq,ir))
            for i in range(1,len(lines)):
                fh.writelines("%s" % lines[i]) #write back
    except IOError:
        print("ERROR adding comment to output file %s" % fn)
        sys.exit(1)


def get_epsilon_optics_vasp(outcarfn,freq):
    delta=0.5
    try:
        outcar_fh = open(outcarfn, 'r')
    except IOError:
        print("ERROR Couldn't open OUTCAR file, skip...")
    while True:
        line=outcar_fh.readline()
        if not line:
            break
        if ('REAL DIELECTRIC FUNCTION' in line):
            break
#            buf=[]
#            print(line)
    for line in outcar_fh:
        if (strType(line.split()[0]) == 'float'):
            f=float(line.split()[0])
            print(f)
            if(abs(f-freq)<delta):
                print(line)

#                        delta=abs(freq-f)
#                        buf=line.split()[1:]
                        #eps=[f,[float(buf[i]) for i in range(len(buf))]]
#                eps=[line.split()[1], line.split()[4], line.split()[6],
#                     line.split()[4], line.split()[2], line.split()[5],
#                     line.split()[6], line.split()[5], line.split()[3]]
#       0   1  2  3  4  5  6
#     Freq XX YY ZZ  XY YZ ZX
                eps=[float(line.split()[i]) for i in (1,4,6,4,2,5,6,5,3)]
                return(np.array(eps))


def get_epsilon_dfpt_vasp(vasprunfn):
    eps = []
    try:
        for event, root in _iterparse(vasprunfn,'calculation'):
            for element in root:
                if element.tag == 'varray':
                    if 'name' in element.attrib:
                        if element.attrib['name'] == 'epsilon':
                            for vectors in element.findall('./v'):
                                    eps.append([float(vec) for vec in vectors.text.split()])
    except:
        print ('Error reading xml file %s', vasprunfnm)
    return(np.array(eps).flatten())

def get_epsilon_dfpt(basedirname,modenum,calculator):
    if(calculator=='vasp'):
        fnm="%s-%03d-1" %(basedirname, (modenum+1))
        fnp="%s-%03d+1" %(basedirname, (modenum+1))
        if os.path.isfile(os.path.join(fnm,'vasprun.xml')) and os.path.isfile(os.path.join(fnp,'vasprun.xml')):
            try:
                epsilon_m=get_epsilon_dfpt_vasp(os.path.join(fnm,'vasprun.xml'))
                epsilon_p=get_epsilon_dfpt_vasp(os.path.join(fnp,'vasprun.xml'))
#                print('Epsilon is: ', epsilon_m,epsilon_p)
            except:
                epsilon_m=np.zeros(9)
                epsilon_p=np.zeros(9)
        else:
            epsilon_m=np.zeros(9)
            epsilon_p=np.zeros(9)
            print('Outputfile %s does not exist!' % os.path.join(fnm,'vasprun.xml'))
    elif(calculator=='castep'):
        fnm="%s-%03d-1" %(basedirname, (modenum+1))
        fnp="%s-%03d+1" %(basedirname, (modenum+1))
        if os.path.isfile(os.path.join(fnm,'shiftcell.castep')) and os.path.isfile(os.path.join(fnp,'shiftcell.castep')):
            try:
                epsilon_m=get_epsilon_dfpt_castep(os.path.join(fnm,'shiftcell.castep'))
                epsilon_p=get_epsilon_dfpt_castep(os.path.join(fnp,'shiftcell.castep'))
            except:
                epsilon_m=np.zeros(9)
                epsilon_p=np.zeros(9)
        else:
            epsilon_m=np.zeros(9)
            epsilon_p=np.zeros(9)
            print('Outputfile %s does not exist!' % os.path.join(fnm,'shiftcell.castep'))
#    print('Epsilon is: ', epsilon_m,epsilon_p)
    elif(calculator=='abinit'):
        fnm="%s-%03d-1" %(basedirname, (modenum+1))
        fnp="%s-%03d+1" %(basedirname, (modenum+1))
        if os.path.isfile(os.path.join(fnm,'shiftcell.abo')) and os.path.isfile(os.path.join(fnp,'shiftcell.abo')):
            try:
                epsilon_m=get_epsilon_abinit(os.path.join(fnm,'shiftcell.abo'))
                epsilon_p=get_epsilon_abinit(os.path.join(fnp,'shiftcell.abo'))
            except:
                epsilon_m=np.zeros(9)
                epsilon_p=np.zeros(9)
        else:
            epsilon_m=np.zeros(9)
            epsilon_p=np.zeros(9)
            print('Outputfile %s does not exist!' % os.path.join(fnm,'shiftcell.abo'))
    return(epsilon_m,epsilon_p)

#epsm = get_epsilon_optic(os.path.join(fnm,'shiftcell_epsilon.dat'))
def get_epsilon_optics(basedirname,modenum,calculator,freq=0.0,cvol=1.0):
    if(calculator=='vasp'):
        fnm="%s-%03d-1" %(basedirname, (modenum+1))
        fnp="%s-%03d+1" %(basedirname, (modenum+1))
        if os.path.isfile(os.path.join(fnm,'OUTCAR')) and os.path.isfile(os.path.join(fnp,'OUTCAR')):
            try:
                epsilon_m=get_epsilon_optics_vasp(os.path.join(fnm,'OUTCAR'),freq)
                epsilon_p=get_epsilon_optics_vasp(os.path.join(fnp,'OUTCAR'),freq)
            except:
                print('Filed to extract epsilon values')
                epsilon_m=np.zeros(9)
                epsilon_p=np.zeros(9)
        else:
            epsilon_m=np.zeros(9)
            epsilon_p=np.zeros(9)
            print('Outputfile %s does not exist!' % os.path.join(fnm,'OUTCAR'))

    elif(calculator=='castep'):
        fnm="%s-%03d-1" %(basedirname, (modenum+1))
        fnp="%s-%03d+1" %(basedirname, (modenum+1))
        if os.path.isfile(os.path.join(fnm,'shiftcell_epsilon.dat')) and os.path.isfile(os.path.join(fnp,'shiftcell_epsilon.dat')):
            try:
                epsilon_m=get_epsilon_optics_castep(os.path.join(fnm,'shiftcell_epsilon.dat'))
                epsilon_p=get_epsilon_optics_castep(os.path.join(fnp,'shiftcell_epsilon.dat'))
            except:
                epsilon_m=np.zeros(9)
                epsilon_p=np.zeros(9)
        else:
            epsilon_m=np.zeros(9)
            epsilon_p=np.zeros(9)
            print('Outputfile %s does not exist!' % os.path.join(fnm,'shiftcell_epsilon.dat'))
    elif(calculator=='cp2k'):
        fnm="%s-%05d-1" %(basedirname, (modenum+1))
        fnp="%s-%05d+1" %(basedirname, (modenum+1))
        if os.path.isfile(os.path.join(fnm,'polar.out')) and os.path.isfile(os.path.join(fnp,'polar.out')):
            try:
                epsilon_m=get_epsilon_cp2k(os.path.join(fnm,'polar.out'))/cvol*Angst2Bohr**3
                epsilon_p=get_epsilon_cp2k(os.path.join(fnp,'polar.out'))/cvol*Angst2Bohr**3
            except:
                epsilon_m=np.zeros(9)
                epsilon_p=np.zeros(9)
        else:
            epsilon_m=np.zeros(9)
            epsilon_p=np.zeros(9)
            print('Outputfile %s does not exist!' % os.path.join(fnm,'polar.out'))
    elif(calculator=='crystal'):
        fnm="%s-%05d-1" %(basedirname, (modenum+1))
        fnp="%s-%05d+1" %(basedirname, (modenum+1))
        if os.path.isfile(os.path.join(fnm,"".join('%s.out' %fnm))) and os.path.isfile(os.path.join(fnp,"".join('%s.out' %fnp))):
            try:
                epsilon_m=get_epsilon_crystal(os.path.join(fnm,"".join('%s.out' %fnm)))/cvol*Angst2Bohr**3
                epsilon_p=get_epsilon_crystal(os.path.join(fnp,"".join('%s.out' %fnp)))/cvol*Angst2Bohr**3
            except:
                epsilon_m=np.zeros(9)
                epsilon_p=np.zeros(9)
        else:
            epsilon_m=np.zeros(9)
            epsilon_p=np.zeros(9)
            print('Outputfile %s does not exist!' % os.path.join(fnm,'polar.out'))
    elif(calculator=='cp2kv6'):
        fnm="%s-%05d-1" %(basedirname, (modenum+1))
        fnp="%s-%05d+1" %(basedirname, (modenum+1))
        if os.path.isfile(os.path.join(fnm,'polar.out')) and os.path.isfile(os.path.join(fnp,'polar.out')):
                epsilon_m=get_epsilon_cp2kv6(os.path.join(fnm,'polar.out'))/cvol*Angst2Bohr**3
                epsilon_p=get_epsilon_cp2kv6(os.path.join(fnp,'polar.out'))/cvol*Angst2Bohr**3
        else:
            epsilon_m=np.zeros(9)
            epsilon_p=np.zeros(9)
            print('Outputfile %s does not exist!' % os.path.join(fnm,'polar.out'))

    return(epsilon_m,epsilon_p)

def genxyz(fn,modenum,freq,basis,species,cart):
    print("Generating file %s" % fn)
    try:
        fh = open(fn, 'w')
    except IOError:
        print("ERROR Couldn't open output file %s for writing" % fn)
        return(-1)
    fh.write('%d\n' % len(species))
    fh.write('%s\n' % 'Epsilon calculation for %d mode, %f cm-1' % (modenum,freq))
    for i in range(len(species)):
        fh.write("%s  %s\n" % (species[i],
                             "".join("    % 15.10f" % c for c in cart[i])))

def gengui(fn,modenum,freq,basis,numbers,cart,ir=''):
    print("Generating file %s" % fn)
    try:
        fh = open(fn, 'w')
    except IOError:
        print("ERROR Couldn't open output file %s for writing" % fn)
        return(-1)
    fh.write(' 3 1 3 %s\n' % 'Epsilon calculation for %d mode, %f cm-1 IR: %s' % (modenum,freq,ir))
    for i in range(3):
        fh.write('%s\n' % "".join("%12.9f " % b for b in basis[i]))
    fh.write(" 1\n 1.0000000 0.0000000 0.0000000\n")
    fh.write(" 0.0000000 1.0000000 0.0000000\n")
    fh.write(" 0.0000000 0.0000000 1.0000000\n")
    fh.write(" 0.0000000 0.0000000 0.0000000\n")
    fh.write(" %d\n" % len(numbers))
    for i in range(len(numbers)):
        fh.write("%d  %s\n" % (numbers[i],
                             "".join("    % 15.10f" % c for c in cart[i])))
def get_epsilon_cp2k(fn):
    epsilon=np.zeros(9)
    try:
        fh = open(fn, 'r')
    except IOError:
        print("ERROR Couldn't open output file %s for reading" % fn)
        return(-1)
    for line in fh:
        if ('Polarizability tensor [a.u.]' in line):
            break
# Components sequence in CP2K output
# xx,yy,zz
# xy,xz,yz
# yx,zx,zy
# Components sequence in epsilon array:
#                     0    1    2    3    4    5    6    7    8
#      9-components:  xx   xy   xz   yx   yy   yz   zx   zy   zz
    line=fh.readline()
    for i in range(3):
        epsilon[i*4]=float(line.split()[i+2])

    line=fh.readline()
    for i in range(2):
        epsilon[1+i]=float(line.split()[i+2])
    epsilon[5]=float(line.split()[4])

    line=fh.readline()
    epsilon[3]=float(line.split()[2])
    for i in range(2):
        epsilon[6+i]=float(line.split()[i+3])   
# Symmitrize
    epsilon[1]=(epsilon[1]+epsilon[3])/2
    epsilon[3]=epsilon[1]
    epsilon[2]=(epsilon[2]+epsilon[6])/2
    epsilon[6]=epsilon[2]
    epsilon[5]=(epsilon[5]+epsilon[7])/2
    epsilon[7]=epsilon[5]
#    print(epsilon)

    return(epsilon)

def get_epsilon_cp2kv6(fn):
    epsilon=np.zeros(9)
    print('cp2kv6')
    try:
        fh = open(fn, 'r')
    except IOError:
        print("ERROR Couldn't open output file %s for reading" % fn)
        return(-1)
    for line in fh:
        if ('POLARIZABILITY TENSOR (atomic units)' in line):
            print(line)
            break

    line=fh.readline()
    for i in range(3):
        epsilon[i*4]=float(line.split()[i+1])

    line=fh.readline()
    for i in range(2):
        epsilon[1+i]=float(line.split()[i+1])
    epsilon[5]=float(line.split()[3])

    line=fh.readline()
    epsilon[3]=float(line.split()[1])
    for i in range(2):
        epsilon[6+i]=float(line.split()[i+2])   
# Symmitrize
    epsilon[1]=(epsilon[1]+epsilon[3])/2
    epsilon[3]=epsilon[1]
    epsilon[2]=(epsilon[2]+epsilon[6])/2
    epsilon[6]=epsilon[2]
    epsilon[5]=(epsilon[5]+epsilon[7])/2
    epsilon[7]=epsilon[5]
#    print(epsilon)

    return(epsilon)

def get_epsilon_crystal(fn):
    epsilon=np.zeros(9)
    try:
        fh = open(fn, 'r')
    except IOError:
        print("ERROR Couldn't open output file %s for reading" % fn)
        return(-1)
    for line in fh:
        if ('SUSCEPTIBILITY (CHI(1)) TENSORS' in line):
            break
    for line in fh:
        if ('ALPHA(REAL,' in line):
            break
# Components sequence in Crystal output
# COMPONENT    ALPHA(REAL, IMAGINARY)         EPSILON       CHI(1) 
#  XX      602.1108        0.0000        7.9693        6.9693    
#  XY       -0.0000        0.0000       -0.0000       -0.0000    
#  XZ       -0.0000        0.0000       -0.0000       -0.0000    
#  YY      564.6560        0.0000        7.5358        6.5358    
#  YZ        0.0000        0.0000        0.0000        0.0000    
#  ZZ      375.7785        0.0000        5.3496        4.3496    
# Components sequence in epsilon array:
#                     0    1    2    3    4    5    6    7    8
#      9-components:  xx   xy   xz   yx   yy   yz   zx   zy   zz

    for i in range(3):
        line=fh.readline()
        epsilon[i]=float(line.split()[1])
    epsilon[3]=epsilon[1]
    epsilon[6]=epsilon[2]
    line=fh.readline()
    epsilon[4]=float(line.split()[1])
    line=fh.readline()
    epsilon[5]=float(line.split()[1])
    epsilon[7]=epsilon[5]
    line=fh.readline()
    epsilon[8]=float(line.split()[1])

#    print(epsilon)

    return(epsilon)
def genabinit(fn,modenum,freq,basis,species,cart,ir=''):
    cell = Atoms(cell=basis, symbols=species, 
            scaled_positions=cart2direct(cart,basis))
    try:
        write_abinit(fn, cell)
    except:
        print('Error writing file %s' % fn)
    try:
        with open(fn,"r") as fh:
            lines = fh.readlines() #read
        with open(fn, "w") as fh:
            fh.write('#%s\n' % 'Epsilon calculation for %d mode, %f cm-1 IR %s' % (modenum,freq,ir))
            for i in range(len(lines)):
                fh.writelines("%s" % lines[i]) #write back
    except IOError:
        print("ERROR adding comment to output file %s" % fn)
        sys.exit(1)

def get_epsilon_abinit(fn):
    e=[]
    try:
        abinit_fh = open(fn, 'r')
    except IOError:
        print ("ERROR Couldn't open abinit output file %s, exiting...\n" % abinitfn)
        sys.exit(1)

    while True:
        line=abinit_fh.readline()
        if not line: break
        if 'Dielectric tensor, in cartesian coordinates' in line:
            while True:
                sline=abinit_fh.readline()
                if not sline: break
                if 'Effective charges' in sline: break
                if (re.match('\s*\d+\s*\d+',sline)):
                    e.append(float(sline.split()[4]))
            break
    return(np.array(e))

Angst2Bohr=1.889725989
#print sqrt(hbar/AMU/10e12)*10e10 #Angstrom

parser = argparse.ArgumentParser(description='The program is to calculate Raman tensor with finite displacements method with VASP/CASTEP/CP2K/ABINIT/CRYSTAL code')

parser.add_argument("-i", "--input", action="store", type=str, dest="str_fn", help="Input filename with structure (POSCAR/castep.cell/input.inp/input.abi)")
parser.add_argument("-o", "--output", action="store", type=str, dest="out_fn", help="Output filename")
#parser.add_argument("-d", "--dynmat", action="store", type=str, dest="dynmat_fn", default='qpoints.yaml', help="Dynmat in yaml format filename")
parser.add_argument("-f", action="store", type=str, dest="fsetfn", default='FORCE_SETS', help="Force_sets filename")
parser.add_argument("--readfc", dest="read_force_constants", action="store_true",
                                                    default=False, help="Read FORCE_CONSTANTS")
parser.add_argument("-s", "--soft", dest="calc", action="store",
                                                    help="Calculator software vasp/castep/cp2k/ABINIT/crystal")
parser.add_argument("-p", "--policy", action="store", type=str, dest="policy", default='displ',
  help="Script modes. 'displ' -- generate input files; 'calcdfpt' -- calculate raman tensor (Epsilon calculated with DFPT method); calc --calculate raman tensor (optics/linear response method)")
parser.add_argument("-D", "--delta", action="store", type=float, dest="delta", default=0.1, help="Shift vector Delta")
parser.add_argument("-m", "--mult", action="store", type=float, dest="mult", default=1.0, help="Intensity multiplier")
parser.add_argument("--freq", action="store", type=float, dest="epsfreq", default=0.0, help="Frequency for epsilon delta. Default 0")
parser.add_argument("--irreps", dest="irreps", action="store_true",
                                                    default=False, help="Find Irreducible representations and print in input files")


args = parser.parse_args()

if (args.str_fn == None):
    print('Error. No input filename was given.')
    sys.exit(1)

basedirname='EPSILON'


if(args.calc==None):
    print('Error. Calculator name is missed')
    sys.exit(1)
calc=""
if(args.calc=="castep"):
    factorcm=521.47083
    calc='castep'
elif(args.calc=="vasp"):
    factorcm=521.47083
    calc='vasp'
elif(args.calc=="cp2k"):
    factorcm=3739.4256800756
    calc='cp2k'
elif(args.calc=="cp2kv6"):
    factorcm=3739.4256800756
    calc='cp2k'
elif(args.calc=="abinit"):
    factorcm=716.85192105135115965589
    calc='abinit'
elif(args.calc=="crystal"):
    factorcm=521.47083
    calc='crystal'
else:
    print('Wrong calculator name %s' % args.calc)
    sys.exit(1)

if (args.read_force_constants):
    print("Read force contants from FORCE_CONSTANTS file")
    ph = phonopy.load(supercell_matrix=[1, 1, 1],
                  primitive_matrix=[1, 0, 0, 0, 1, 0, 0, 0, 1],
                  unitcell_filename=args.str_fn,
                  calculator=calc, factor=factorcm,
                  force_constants_filename='FORCE_CONSTANTS')
else:
    ph = phonopy.load(supercell_matrix=[1, 1, 1],
                  primitive_matrix=[1, 0, 0, 0, 1, 0, 0, 0, 1],
                  unitcell_filename=args.str_fn,
                  calculator=calc, factor=factorcm,
                  force_sets_filename=args.fsetfn)

species=ph.primitive.get_chemical_symbols()
numbers=ph.primitive.get_atomic_numbers()
masses=ph.primitive.get_masses()
natom=ph.primitive.get_number_of_atoms()
basis=ph.primitive.get_cell()
if(args.calc=='castep'):
    magmoms=ph.primitive.get_magnetic_moments()
print("primitive")
print(ph.primitive.get_scaled_positions())
cart=direct2cart(ph.primitive.get_scaled_positions(),basis)
cvol=np.dot(basis[0],np.cross(basis[1],basis[2]))
print(species)
print(numbers)
print(cart)
print(cart2direct(cart,basis))
print(basis)

if (args.irreps):
# Hack to do irrep analysis. Spin order will be restored later
    ph.primitive.set_magnetic_moments(None)
# Set IR for Gamma point
    ph.set_irreps([0.0,0.0,0.0])
    ir_labels=[]
    print(ph.get_irreps()._ir_labels)
    for ir in ph.get_irreps()._ir_labels:
        if ( ir is None):
            ir_labels.append('None')
        elif ('T' in ir):
            for i in range(3):
                ir_labels.append(ir)
            continue
        elif  ('E' in ir):
            for i in range(2):
                ir_labels.append(ir)
            continue
        else:
            ir_labels.append(ir)
    print(ir_labels)
else:
    ir_labels=["" for x in range(natom*3)]

dmat = ph.get_dynamical_matrix_at_q([0,0,0])
eigvals, evecs = np.linalg.eigh(dmat)

frequencies=np.sqrt(np.abs(eigvals.real)) * np.sign(eigvals) 
for i in range(len(frequencies)):
# Accoustic modes could be zero, shif frequency a little bit
    if(abs(frequencies[i]*factorcm)<1):
        frequencies[i]=1/factorcm

eigvecs=evecs.T

#for i in range(len(frequencies)):
#    print(frequencies[i]*factorcm)
#    for j in range(len(species)):
#        print("".join("% 9.7f " % (eigvecs[i][j*3+k]/sqrt(masses[j])) for k in range(3)))

# GENERATION OF DISPLACEMENTS INPUT FILES
#j - mode number; i -atom number

if (args.policy == 'displ'):
    for i in range(len(frequencies)):
        cartshiftdm=[]
        cartshiftdp=[]

#        print ('mode: %d %8.5f' % ((i+1),frequencies[i]))

#        print('====eigenvector:====')
#        for j in range(len(species)):
#            print(" ".join('% 9.7f' % eigvecs[i][j*3+l].real for l in range(3)))
#        print('delta = %f' % args.delta)

#        print('shiftvector:')

        for j in range(len(species)):
            shiftvec=[0.0e0,0.0e0,0.0e0]
            for l in range(3):
                shiftvec[l]=eigvecs[i][j*3+l].real*args.delta*18.362*sqrt(1/(masses[j]*abs(frequencies[i])*factorcm))/Angst2Bohr
#            print(["%10.7f" % shiftvec[l] for l in range(3)])
            cartshiftdm.append(cart[j]-np.array(shiftvec))
            cartshiftdp.append(cart[j]+np.array(shiftvec))

#        print('shifted:')
#        for cartat in cartshiftdm:
#            print ("%12.9f %12.9f %12.9f" % (cartat.tolist()[0], cartat.tolist()[1], cartat.tolist()[2]))
        if(args.calc=='castep'):
            fnm="shiftcell-%03d-1.cell" % (i+1)
            fnp="shiftcell-%03d+1.cell" % (i+1)
            gencastep(fnm,species,basis,cartshiftdm, magmoms=magmoms, 
                    commentstr='freq = %9.4f cm-1; Delta=%6.4f' % ((frequencies[i] * factorcm),args.delta),ir=ir_labels[i])
            gencastep(fnp,species,basis,cartshiftdp, magmoms=magmoms, 
                    commentstr='freq = %9.4f cm-1; Delta=%6.4f' % ((frequencies[i] * factorcm),args.delta),ir=ir_labels[i])
        elif(args.calc=='vasp'):
            poscarfnm="POSCAR-%03d-1" % (i+1)
            poscarfnp="POSCAR-%03d+1" % (i+1)
            genposcar(poscarfnm,i+1,frequencies[i]*factorcm,basis,species,cartshiftdm,ir=ir_labels[i])
            genposcar(poscarfnp,i+1,frequencies[i]*factorcm,basis,species,cartshiftdp,ir=ir_labels[i])
#Hack to support cp2k v6
        elif(calc=='cp2k'):
            poscarfnm="shiftcell-%05d-1.xyz" % (i+1)
            poscarfnp="shiftcell-%05d+1.xyz" % (i+1)
            genxyz(poscarfnm,i+1,frequencies[i]*factorcm,basis,species,cartshiftdm)
            genxyz(poscarfnp,i+1,frequencies[i]*factorcm,basis,species,cartshiftdp)
        elif(args.calc=='abinit'):
            fnm="shiftcell-%03d-1.abi" % (i+1)
            fnp="shiftcell-%03d+1.abi" % (i+1)
            genabinit(fnm,i+1,frequencies[i]*factorcm,basis,species,cartshiftdm,ir=ir_labels[i])
            genabinit(fnp,i+1,frequencies[i]*factorcm,basis,species,cartshiftdp,ir=ir_labels[i])
        elif(calc=='crystal'):
            fnm="shiftcell-%05d-1.gui" % (i+1)
            fnp="shiftcell-%05d+1.gui" % (i+1)
            gengui(fnm,i+1,frequencies[i]*factorcm,basis,numbers,cartshiftdm,ir=ir_labels[i])
            gengui(fnp,i+1,frequencies[i]*factorcm,basis,numbers,cartshiftdp,ir=ir_labels[i])
    if(args.calc=='cp2k'):
        print('A B C for input file is:')
        print('A %s' % "".join("%14.9f" % b for b in basis[0]))
        print('B %s' % "".join("%14.9f" % b for b in basis[1]))
        print('C %s' % "".join("%14.9f" % b for b in basis[2]))

######################## Calculate RAMAN ###########################
else:
    if (args.out_fn == None):
        print('Error. No output filename was given.')
        sys.exit(1)

    try:
        out_fh = open(args.out_fn, 'w')
    except IOError:
        print("ERROR Couldn't open output file for writing, exiting...")
        sys.exit(1)
    out_fh.write("# N     Freq         xx         xy         xz         yx         yy         yz         zx         zy         zz         Alpha         Gamma2         Ipar         Iperp         Itot\n")
    for i in range(len(frequencies)):
        cartshiftdm=[]
        cartshiftdp=[]
        G0=0
        G1=0
#        print('mode number: %d; dirname %s; freq: %f' %((j+1),vasprunfnp,frequencies[natom*3-j-1]*args.freqfactor))
########################  DFPT  ###########################
        if (args.policy == 'calcdfpt'):
            epsm, epsp = get_epsilon_dfpt(basedirname,i,args.calc)
########################  OPTIC  ###########################
        else:
            epsm, epsp = get_epsilon_optics(basedirname,i,args.calc,freq=args.epsfreq,cvol=cvol)

        if (np.linalg.norm(epsp)==0 or np.linalg.norm(epsm)==0):
            print ('Epsilon data damaged in for mode %d. Norm is zero' % (i+1))
            continue
        print('Got epsilon values difference:')

# Atomic units: sqrt(Bohr/amu) (me=1, e=1,hbar=1)
        if (args.calc=='abinit'):
            alpha=(epsp-epsm)*sqrt(cvol)/(4*pi)/(2*args.delta*18.362*sqrt(1/(abs(frequencies[i])*factorcm)))*sqrt(Angst2Bohr**2*9.10938356/1.6605402*10**-4)*args.mult
        else:
            alpha=(epsp-epsm)*sqrt(cvol)/(4*pi)/(2*args.delta*18.362*sqrt(1/(abs(frequencies[i])*factorcm)))*sqrt(Angst2Bohr*9.10938356/1.6605402*10**-4)*args.mult
# Calculate invariance in Long's notation
#                     0    1    2    3    4    5    6    7    8
#      9-components:  xx   xy   xz   yx   yy   yz   zx   zy   zz
        Alpha=((alpha[0]+alpha[4]+alpha[8]))/3
        Gamma2=( (alpha[0]-alpha[4])**2+(alpha[4]-alpha[8])**2+(alpha[8]-alpha[0])**2 )/2
        Gamma2+=3*( (alpha[1])**2 + (alpha[2])**2 + (alpha[5])**2 )

        Iperp=Gamma2/15
        Ipar=( 45*Alpha**2 + 4*Gamma2 )/45
        print('%4d\t%8.2f\t%12.8f\t%12.8f\t%12.8f' % ((i+1), frequencies[i]*factorcm,Ipar,Iperp,(Ipar+Iperp)))

        out_fh.write('%4d  %8.2f  ' % ((i+1), frequencies[i]*factorcm))
        out_fh.write(' '.join(" %9.6f" % alpha[i] for i in range(9)))
        out_fh.write('  % 12.10f % 12.10f % 12.10f % 12.10f % 12.10f\n' % (Alpha**2,Gamma2,Ipar,Iperp,(Ipar+Iperp)))

    #        print ('\nItotal=%5.3f Iparal=%5.3f Iperp=%5.3f' % ((10*G0+7*G2+5*G1), (10*G0+4*G2), (5*G1+3*G2)))
    # Multiplyer for Intensity at room temperature with 514.5nm excitation line
    #        print ('Itot*C=%10.8f\n' % ( ( (19436.35-frequencies[j])**4 )  / ( 1 - exp ( -0.2281*Temp*frequencies[j] ) ) / (30*1E12*frequencies[j])*(10*G0+7*G2+5*G1)   )   )


    out_fh.close()
    
