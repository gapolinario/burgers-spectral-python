#!/bin/bash

# runs Burgers equation with zero forcing
# and a sine initial condition
# then runs the same code with PyPDE and
# and produce some figures to compare the output
# result is deterministic

# from main folder, run
# bash tests/test_burgers_initial_vs_pypde.sh
# on btphm4 takes 20s
# then go to Test Burgers PyPDE.ipynb and run it

# N Ltotal Lrelative Ttotal nu sqeps NTsave CFL_const
params='8 1. 3  1. 2 1. 100  0.1'
# nlinear fkernel initial_value scheme saveformat
textparams='True zero_forcing sine_initial ETD real'
prefix=90000

(set -x; # prints to the screen the command being run
python3 -O scripts/burgers.py $prefix $params $textparams
)
