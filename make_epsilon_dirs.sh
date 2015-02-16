#!/bin/bash

INITDIR=$1
TAIL=$2
CWD=`pwd`

if [[ ! -d $INITDIR ]]; then
    echo "Error init directory with abinit input files does not exist."
    echo "Usage $0 init_dir abinit_tail.in"
    exit 1
fi

if [[ ! -f $TAIL ]]; then
    echo "Error abinit_Tail file does not exist."
    echo "Usage $0 init_dir abinit_tail.in"
    exit 1
fi

if [[ ! -f $INITDIR/shiftcell.files ]]; then
    echo "Error. You should prepare init directory carefully. Make shure shiftcell.files file is exist"
    exit 1
fi

for file in $INITDIR/shiftcell-*; do
	fn=${file%.in}
	num=${fn#*-}
	if [[ -d "epsilon-$num" ]]; then
	    echo "Error directory epsilon-$num already exist"
	    exit 1
	fi
	echo "Making epsilon-$num directory"
        mkdir epsilon-$num
	cp $INITDIR/*.fhi epsilon-$num/
	cp $INITDIR/*.xml epsilon-$num/
	cp $INITDIR/shiftcell.files epsilon-$num/
	cp $file epsilon-$num/shiftcell.in

	echo "Generating new input using ${TAIL} pattern"

	echo -e "\n" >> "epsilon-$num"/shiftcell.in
	cat $TAIL >> "epsilon-$num"/shiftcell.in
done
