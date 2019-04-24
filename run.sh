#!/bin/bash
export CONFIG="$1"
export CHARGE="$2"
jdftx -i totalE.in | tee $CONFIG-totalE.out
jdftx -i vibrations.in | tee $CONFIG-vibrations.out
