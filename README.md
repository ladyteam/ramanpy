# ramanpy
Raman tensor calculation using finite displacement method

Python script to calculate Raman tensor by finite displacement method. The script based on PHONOPY project.
http://phonopy.sourceforge.net/ by Atsushi Togo.

The goal is to calculate Raman tensor with ABINIT package (http://www.abinit.org). The idea is to calculate 
static dielectric tensor for shifted cells. The procedure is as follows:

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

9. Now you have a set of directories with name pattern epsilon-NNN, there NNN -- is a mode number.
   Run calculations in appropriative directories (usially one wants to skip first 3 acoustic modes)

10. When calculations is finished run the script with the following command to plot raman tensors:
  ./raman_abinit.py -i PHONON/supercell.in -p calc
