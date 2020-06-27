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
# DOI: 
#
# The code is based on PHONOPY API https://phonopy.github.io/phonopy/
#
# Author: Eugene Roginskii

import h5py
import matplotlib.pyplot as plt 
import phonopy
from phonopy.interface.castep import write_castep
from phonopy.structure.atoms import PhonopyAtoms as Atoms
#from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
#from phonopy.phonon.band_structure import get_band_structure_dict
import argparse
import matplotlib
import numpy as np
import sys
from math import sqrt
import re
from math import (
        sqrt, pi, exp)
from shutil import move
import os
import datetime
import time



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
        direct.append(np.dot(atcart,np.linalg.inv(np.transpose(basis))))
    return np.array(direct)


def get_epsilon_dfpt(castep_fn):
    e=[]
    try:
        fh = open(castep_fn, 'r')
    except IOError:
        print ("ERROR Couldn't open abinit output file %s, exiting...\n" % abinitfn)
        sys.exit(1)

    while True:
        line=fh.readline()
        if not line: break
        if 'Dielectric tensor, in cartesian coordinates' in line:
            while True:
                sline=abinit_fh.readline()
                if not sline: break
                if 'Effective charges' in sline: break
                if (re.match('\s*\d+\s*\d+',sline)):
                    e.append(float(sline.split()[4]))
            break
    eps=np.array(e).reshape(3,3)
    return eps

def gencastep(fn,species,basis,cart,commentstr):
    cell = Atoms(cell=basis, symbols=species, 
            scaled_positions=cart2direct(cart,basis))
    try:
        write_castep(fn, cell)
    except:
        print('Error writing file %s' % fn)
    try:
        with open(fn,"r") as fh:
            lines = fh.readlines() #read
        with open(fn, "w") as fh:
            fh.write('#%s\n' % commentstr)
            fh.writelines(lines) #write back
    except IOError:
        print("ERROR adding comment to output file %s" % fn)
        sys.exit(1)

####### Extract permittivity from DFPT #######
def get_epsilon_dfpt(fn):
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
def get_epsilon_optic(fn):
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

    return(np.array([e[0],e[3],e[4],e[3],e[1],e[5],e[4],e[5],e[2]]))

# CASTEP freq to cm-1 factor
factorcm=521.47083

Angst2Bohr=1.889725989
basedirname='epsilon'


parser = argparse.ArgumentParser(description='The program is to calculate Raman tensor with finite displacements method with CASTEP code')

parser.add_argument("-i", "--input", action="store", type=str, dest="castep_fn", help="CASTEP input filename (.cell)")
parser.add_argument("-o", "--output", action="store", type=str, dest="out_fn", help="Output filename")
parser.add_argument("-f", action="store", type=str, dest="fsetfn", default='FORCE_SETS', help="Force_sets filename")
parser.add_argument("-p", "--policy", action="store", type=str, dest="policy", default='displ',
  help="Script modes. 'displ' -- generate input files; 'calcdfpt' -- calculate raman tensor (Epsilon calculated with DFPT method); calc --calculate raman tensor (optics method)")
parser.add_argument("-D", "--delta", action="store", type=float, dest="delta", default=0.1, help="Shift vector Delta")
parser.add_argument("-m", "--mult", action="store", type=float, dest="mult", default=1.0, help="Intensity multiplier")

args = parser.parse_args()

if (args.castep_fn == None):
    print('Error. No input filename was given.')
    sys.exit(1)
if (args.out_fn == None):
    print('Error. No output filename was given.')
    sys.exit(1)


ph = phonopy.load(supercell_matrix=[1, 1, 1],
                  primitive_matrix=[1, 0, 0, 0, 1, 0, 0, 0, 1],
                  unitcell_filename=args.castep_fn,
                  calculator='castep', factor=factorcm,
                  force_sets_filename=args.fsetfn)

species=ph.primitive.get_chemical_symbols()
masses=ph.primitive.get_masses()
natom=ph.primitive.get_number_of_atoms()
basis=ph.primitive.get_cell()
cart=direct2cart(ph.primitive.get_scaled_positions(),basis)
cvol=np.dot(basis[0],np.cross(basis[1],basis[2]))
print(species)
print(cart)
print(cart2direct(cart,basis))
print(basis)

dmat = ph.get_dynamical_matrix_at_q([0,0,0])
eigvals, evecs = np.linalg.eigh(dmat)

frequencies=np.sqrt(np.abs(eigvals.real)) * np.sign(eigvals) 
for i in range(len(frequencies)):
# Accoustic modes could be zero, shif frequency a little bit to exculde zero division
    if(abs(frequencies[i]*factorcm)<1):
        frequencies[i]=1/factorcm

eigvecs=evecs.T

if (args.policy == 'displ'):
    for i in range(len(frequencies)):
        cartshiftdm=[]
        cartshiftdp=[]

        print ('mode: %d %8.5f' % ((i+1),frequencies[i]))

        print('====eigenvector:====')
        for j in range(len(species)):
            print(" ".join('% 9.7f' % eigvecs[i][j*3+l].real for l in range(3)))
        print('delta = %f' % args.delta)

        print('shiftvector:')

        for j in range(len(species)):
            shiftvec=[0.0e0,0.0e0,0.0e0]
            for l in range(3):
                shiftvec[l]=eigvecs[i][j*3+l].real*args.delta*18.362*sqrt(1/(masses[j]*abs(frequencies[i])*factorcm)) 
            print(["%10.7f" % shiftvec[l] for l in range(3)])
            cartshiftdm.append(cart[j]-np.array(shiftvec))
            cartshiftdp.append(cart[j]+np.array(shiftvec))

        print('shifted:')
        for cartat in cartshiftdm:
            print ("%12.9f %12.9f %12.9f" % (cartat.tolist()[0], cartat.tolist()[1], cartat.tolist()[2]))
        fnm="shiftcell-%03d-1.cell" % (i+1)
        fnp="shiftcell-%03d+1.cell" % (i+1)
        gencastep(fnm,species,basis,cartshiftdm,'freq = %9.4f cm-1; Delta=%6.4f' % ((frequencies[i] * factorcm),args.delta))
        gencastep(fnp,species,basis,cartshiftdp,'freq = %9.4f cm-1; Delta=%6.4f' % ((frequencies[i] * factorcm),args.delta))

######################## Calculate RAMAN ###########################

else:
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
        fnm="%s-%03d-1" %(basedirname, (i+1))
        fnp="%s-%03d+1" %(basedirname, (i+1))
        if os.path.isfile(os.path.join(fnm,'shiftcell.castep')) and os.path.isfile(os.path.join(fnp,'shiftcell.castep')):
            try:
########################  DFPT  ###########################
                if (args.policy == 'calcdfpt'):
                    epsm = get_epsilon_dfpt(os.path.join(fnm,'shiftcell.castep'))
                    epsp = get_epsilon_dfpt(os.path.join(fnp,'shiftcell.castep'))
########################  OPTIC  ###########################
                else:
                    if os.path.isfile(os.path.join(fnm,'shiftcell_epsilon.dat')) and os.path.isfile(os.path.join(fnp,'shiftcell_epsilon.dat')):
                        epsm = get_epsilon_optic(os.path.join(fnm,'shiftcell_epsilon.dat'))
                        epsp = get_epsilon_optic(os.path.join(fnp,'shiftcell_epsilon.dat'))
                    else:
                        print('Error. No output optic files exist.')
                        sys.exit(1)

                if (len(epsm)<9):
                    print ('epsilon data damaged in file %s' % fnm)
                    continue
                if (len(epsp)<9):
                    print ('epsilon data damaged in file %s' % fnp)
                    continue

# Atomic units: sqrt(Bohr/amu) (me=1, e=1,hbar=1)
                alpha=(epsm-epsp)*sqrt(cvol)/(4*pi)/(2*args.delta*18.362*sqrt(1/(abs(frequencies[i])*factorcm)))*sqrt(Angst2Bohr*9.10938356/1.6605402*10**-4)*args.mult

# Calculate invariance in Long's notation
                Alpha=((alpha[0]+alpha[4]+alpha[8]))/3
                Gamma2=( (alpha[0]-alpha[4])**2+(alpha[4]-alpha[8])**2+(alpha[8]-alpha[0])**2 )/2
                Gamma2+=3*( (alpha[1])**2 + (alpha[2])**2 + (alpha[5])**2 )

                Iperp=Gamma2/15
                Ipar=( 45*Alpha**2 + 4*Gamma2 )/45
                print('%4d\t%8.2f\t%12.8f\t%12.8f\t%12.8f' % ((i+1), frequencies[i]*factorcm,Ipar,Iperp,(Ipar+Iperp)))

                out_fh.write('%4d  %8.2f  ' % ((i+1), frequencies[i]*factorcm))
                out_fh.write(' '.join(" %9.6f" % alpha[i] for i in range(9)))
                out_fh.write('  % 12.10f % 12.10f % 12.10f % 12.10f % 12.10f\n' % (Alpha**2,Gamma2,Ipar,Iperp,(Ipar+Iperp)))
            except:
                print('Filed to calculate raman for mode %d' % (i+1))


        else:
            print("No data for mode %d in directory %s or %s" % ((i+1),fnm,fnp))

    out_fh.close()
