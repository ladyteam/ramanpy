#!/bin/bash 

mkdir notactive

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
