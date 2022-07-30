#!/bin/bash

TAIL=$1
CWD=`pwd`


for file in shiftcell-*; do
	fn=${file%.in}
	num=${fn#*-}
	if [[ -d "epsilon-$num" ]]; then
	    echo "Error directory epsilon-$num already exist"
	    exit 1
	fi
	echo "Making epsilon-$num directory"
        mkdir EPSILON-$num
	cp $file EPSILON-$num/shiftcell.abi

	echo "Generating new input using ${TAIL} pattern"

	echo -e "\n" >> "EPSILON-$num"/shiftcell.abi
	cat $TAIL >> "EPSILON-$num"/shiftcell.abi
done
