# ramanpy
Raman tensor calculation using finite displacement method

Python script to calculate Raman tensor by finite displacement method. The script based on PHONOPY project.
http://phonopy.sourceforge.net/ by Atsushi Togo.

Warning! The general script is "raman_fd.py" others are obsolete and no more supported.

Instructions (VASP calculator is used in this example. Calculator is given by passing it's name with --soft parameter):
1. Make PHONON directory and calculate Dynamical matrix (IBRION=6 or IBRION=8) then construct FORCE_CONSTANTS file with the command:

$  phonopy --fc vasprun.xml

then go to higher level direcory
cd ..

2. mkdir RAMAN; cd RAMAN; ln -s ../PHONON

3. Gernerate POSCARs files with following command.
$ ./raman_fd.py -i PHONON/POSCAR --readfc --soft vasp -p displ --irreps

help on parameters can be obtained with -h or --help key

4. Edit INCAR file to check the values of variables (NBANDS parameter depends on the system. Normally 
it's something like number of valence band + 250):

PREC = Accurate
LOPTICS = .TRUE.
NBANDS=432
NEDOS=4000
LREAL = .FALSE.
LWAVE = .FALSE.
LCHARG = .FALSE.
IBRION=-1
NSW = 0

5. If irreps option works properly each generated POSCAR will contain Irreducible representations in 
first string as a comment. So check each generates POSCAR-s and delete the ones which are not Raman 
active. For example Au modes are Raman silent, therefore no need to calculate Raman tensor for these modes.* See below
Then make directories for calculations by running bash script make_epsilon_dirs.sh 
(chmod +x make_epsilon_dirs.sh before, to make it executable)

6. Run calculations

7. When finished clear the stuff:

for i in EPSILON-*; do cd $i; rm WAVEDER DOSCAR slurm*out; cd ..; done

8. Now collect the data and save the Raman tensor to file:

$ ./raman_fd.py --readfc -i PHONON/POSCAR -p calc -o raman.dat --soft vasp

Voila check raman.dat for content.


* In order to automitize selecting of Raman active modes simple bash script can be used to filter Raman nonactive modes:
#!/bin/bash 

mkdir notactive

# Set list of nonactive modes for your structure by editing variable "a" below
a=(Au Bu) 
for fn in POSCAR-*; do 
    notactive="False" 
    ir=`head -n1 $fn | awk -F"IR:" '{print($2)}' | xargs`
    for i in ${a[@]}; do 
        if [[ $i == $ir ]]; then 
            echo "Not active mode $fn" 
            mv $fn ./notactive
        fi 
    done
done



Instructions for old scripts:

The goal is to calculate Raman tensor with ABINIT package (http://www.abinit.org) or VASP (https://www.vasp.at/)
by finite difference technique. The idea is to calculate static dielectric tensor for shifted cells. 
The procedure is as follows:

I ABINIT PACKAGE
1. From current directory make PHONON subdir change directory and Calculate force dynamic matrix for 
the structure with PHONOPY.

2. Generate qpoints.yaml file with dynamic matrix at Gamma point. Example of phonopy command:
phonopy --abinit -c supercell.in --qpoints="0 0 0" --factor=716.851928469 --dim="1 1 1" --writedm

3. Go to the upperlevel directory (cd ..)

4. link qpoints.yaml file (ln -s PHONOPY/qpoints.yaml qpoints.yaml)

5. Generate shifted input files by running the script:
  ./raman_abinit.py -i PHONON/supercell.in -p displ

6. Prepare tail.in file (the one consists of task sets). See example tail.in.

7. prepare 'shiftcell' directory and put shiftcell.files (the ABINIT .files file with filenames
   of .in file, out files and potential/PAW-datasets filenames). Copy potentail/PAW-datasets and 
   produced shiftcell-NNN.in files in 'shiftcell' directory. Copy 

8. run make_epsilon_dirs.sh with paramiters:
   ./make_epsilon_dirs.sh shiftcell tail.in

9. Now you have a set of directories with name pattern epsilon-NNN+-1, there NNN -- is a mode number.
   Run calculations in appropriative directories (usially one wants to skip first 3 acoustic modes)

10. When calculations is finished run the script with the following command to plot raman tensor table:
  ./raman_abinit.py -i PHONON/supercell.in -d qpoints.yaml -p calc

II VASP PACKAGE
1. From current directory make PHONON subdir change directory and Calculate force dynamic matrix for 
the structure with PHONOPY.

2. Generate qpoints.yaml file with dynamic matrix at Gamma point. Example of phonopy command:
phonopy  -c POSCAR --qpoints="0 0 0" --factor=521.47083 --dim="1 1 1" --writedm

3. Go to the upperlevel directory (cd ..)

4. link qpoints.yaml file (ln -s PHONOPY/qpoints.yaml qpoints.yaml)

5. Generate shifted input files by running the script:
  ./raman_vasp.py -v PHONON/vasprun.xml -p displ

6. Prepare INPUT file (the one consists of task sets). See example INPUT.epsilon; INPUT.dfpt. And POTCAT and KPOTINS file

7. run make_epsilon_dirs.sh with paramiters:
   ./make_epsilon_dirs_vasp.sh

8. Now you have a set of directories with name pattern EPSILON-NNN+-1, there NNN -- is a mode number.
   Run calculations in appropriative directories (usially one wants to skip first 3 acoustic modes)

9. When calculations is finished run the script with the following command to plot raman tensor table:
  ./raman_vasp.py -v PHONON/vasprun.xml -d qpoints.yaml -p calc
  or 
  ./raman_vasp.py -v PHONON/vasprun.xml -d qpoints.yaml -p calcdfpt
  (depends on the selected method of dielectric response calculation)