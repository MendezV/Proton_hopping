#!/bin/bash


for n_bound in {-2..2}; do
    export nbound="$(echo $n_bound | awk '{printf("%d", 5)}')"
    export ibound="$(echo $n_bound | awk '{printf("%.4f", 0.719643+0.719643*$1/20)}')"
    echo “Outside loop: ${ibound}”
    for (( n_pos=0; n_pos<=${nbound}; n_pos++ )); do
        export inpos="$(echo $n_pos ${ibound} | awk '{printf("%.4f", $2*$1/5  )}')"

	echo ${inpos} 
	#jdftx -i totalE.in
    done
done
