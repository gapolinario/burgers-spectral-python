#!/bin/bash

# runs Burgers equation with stochastic forcing on
# implemented with PyTorch instead of Numpy

# from main folder, run
# bash tests/test_burgers_torch.sh
# can be run in parallel, Fourier transforms benefit from that

# N Ltotal Lrelative Ttotal nu sqeps NTsave cfl_const
params='8 1. 3 10. 2 1. 1000 1.0'
# nlinear fkernel initial_value scheme
textparams='True white_fourier_forcing zero_initial ETD fourier'
prefix=0
ensemble_size=1
numcores=6

ensemble_init=$(($prefix*1000))
ensemble_end=$(($prefix*1000 + $ensemble_size - 1))

(set -x; # prints to the screen the command being run
python3 -O scripts/burgers_tc.py $prefix $params $textparams
)
