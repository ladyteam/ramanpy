#!/bin/bash

CWD=`pwd`

if [[ ! -f INCAR ]]; then
    echo "Error INCAR pattern file does not exist"
    exit 1
fi

if [[ ! -f POTCAR ]]; then
    echo "Error POTCAR file does not exist."
    exit 1
fi

if [[ ! -f KPOINTS ]]; then
    echo "Error KPOINTS file does not exist."
    exit 1
fi


for file in POSCAR-*; do
	num=${file#*-}
	if [[ -d "EPSILON-$num" ]]; then
	    echo "Error directory EPSILON-$num already exist"
	    exit 1
	fi
	echo "Making EPSILON-$num directory"
        mkdir EPSILON-$num
	cp INCAR POTCAR KPOINTS EPSILON-$num/
	cp $file EPSILON-$num/POSCAR
done
