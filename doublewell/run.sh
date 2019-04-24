#!/bin/bash
for n_bound in {4..8}; do
	export n_bound
	export ibound="$(echo $n_bound | awk '{printf("%.4f", 0.719643*$1)}')"
	for n_pos in {1..${n_bound}}; do
    	export inpos="$(echo $n_pos ${ibound} ${n_bound}| awk '{printf("%.4f", 2*$2*$1 /(2*$3) -$2/($3))}')"
    	jdftx -i totalE.in
	done
done
