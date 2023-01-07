This is a pseudospectral implementation of the Burgers equation
3/2 dealiasing is used

# TESTS

1. Deterministic Burgers equation with pseudospectral method
matches the solution obtained with PyPDE

From the main folder, run
`bash tests/test_burgers_initial_vs_pypde.sh`
Then run the jupyter notebook
`tests/Test Burgers Initial PyPDE.ipynb`

2. Ensemble of simulations of the stochastic Heat equation
match analytical solution of the discrete system (small number of Fourier modes)

From the main folder, run
`bash tests/test_discrete_heat_forced.sh`
Then run the jupyter notebook
`tests/Test Discrete Heat Forced.ipynb`

# OUTLINE

1. Deterministic Burgers equation vs PyPDE
2. Stochastic Heat equation ensemble vs discrete solution
3. Stochastic Heat equation ensemble vs continuous solution
4. Stochastic Burgers equation ensemble vs PyPDE
5. Stochastic Burgers equation stationary state
6. Stochastic Burgers equation with Brehier noise
7. Stochastic Burgers equation with JKW algorithm
8. Database for simulations, hashlib indexing
9. Comparison against the Chabanol-Duchon statistics for forced inviscid Burgers turbulence
