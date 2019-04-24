#!/bin/bash


for n_bound in {4..8}; do
    export nbound="$(echo $n_bound | awk '{printf("%d", $1-1)}')"
    export ibound="$(echo $n_bound | awk '{printf("%.4f", 0.719643*$1)}')"
    echo “Outside loop: ${nbound}”
    for (( n_pos=1; n_pos<=${nbound}; n_pos++ )); do
        export inpos="$(echo $n_pos ${ibound} ${nbound}| awk '{printf("%.4f", 2*$2*$1/(2*$3) -$2/($3))}')"

         jdftx -i totalE.in
    done
done
