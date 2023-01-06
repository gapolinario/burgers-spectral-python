#!/bin/bash

# runs Heat equation with stochastic forcing on
# low number of modes and answer is compared against
# analytical solution for discrete problem
# each Fourier mode is independent, since problem is linear
# runs a small ensemble

# dt parameter is 0.1*dx**2
# from main folder, run
# bash tests/test_discrete_heat_forced.sh
# on btphm4 takes 200s with 6 processors
# then go to Test Discrete Heat.ipynb and run it

# N Ltotal Lrelative Ttotal nu sqeps NTsave
params='4 1. 3 10. 2 1. 1000 0.1'
# nlinear fkernel initial_value scheme
textparams='False zero_smooth_fourier_forcing zero_initial ETD fourier'
prefix=91
ensemble_size=300
numcores=6

ensemble_init=$(($prefix*1000))
ensemble_end=$(($prefix*1000 + $ensemble_size - 1))

(set -x; # prints to the screen the command being run
seq $ensemble_init $ensemble_end | xargs -t -I {} -P $numcores python3 scripts/burgers.py {} $params $textparams
)
