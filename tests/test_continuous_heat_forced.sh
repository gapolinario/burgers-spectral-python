#!/bin/bash

# runs Heat equation with stochastic forcing on
# high number of modes and answer is compared against
# analytical solution for continuous problem
# each Fourier mode is independent, since problem is linear
# runs a small ensemble

# from main folder, run
# bash tests/test_continuous_heat_forced.sh
# on btphm4 takes ?? with 6 processors
# then go to Test Continuous Heat.ipynb and run it

# N Ltotal Lrelative Ttotal nu sqeps NTsave
params='9 1. 3 10. 2 1. 1000 1.0'
# nlinear fkernel initial_value scheme
textparams='False zero_smooth_fourier_forcing zero_initial ETD fourier'
prefix=92
ensemble_size=100
numcores=6

ensemble_init=$(($prefix*1000))
ensemble_end=$(($prefix*1000 + $ensemble_size - 1))

(set -x; # prints to the screen the command being run
seq $ensemble_init $ensemble_end | xargs -t -I {} -P $numcores python3 scripts/burgers.py {} $params $textparams
)
