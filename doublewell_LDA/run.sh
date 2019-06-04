#!/bin/bash


for n_bound in {4..8}; do
    export nbound="$(echo $n_bound | awk '{printf("%d", 2*$1)}')"
    export ibound="$(echo $n_bound | awk '{printf("%.4f", 0.719643*$1)}')"
    echo “Outside loop: ${nbound}”
    for (( n_pos=1; n_pos<=${nbound}; n_pos++ )); do
        export inpos="$(echo $n_pos ${ibound} | awk '{printf("%.4f", (2*0.719643)*$1/(2*2) -$2)}')"

	#echo ${inpos} 
	jdftx -i totalE.in
    done
done
