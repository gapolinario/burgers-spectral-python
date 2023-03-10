import torch as tc
from torch.fft import fft,ifft,fftfreq,fftshift
from math import pi
from math import sqrt as m_sqrt
import sys
import os
from json import dump,dumps
from hashlib import sha224
from h5py import File
from numpy import complex128

from _functions import *




######## GET EXTERNAL PARAMETERS ############
#############################################
#############################################
#############################################






# identifier of this realization, external argument
R = int(sys.argv[1])
# spatial parameters
BN = int(sys.argv[2])
Ltotal = float(sys.argv[3])
Lrelative = 1./2**int(sys.argv[4])
Ttotal = float(sys.argv[5]) # total time
nu = 10.**(-int(sys.argv[6]))
sqeps = float(sys.argv[7]) # forcing parameter
NTsave = int(sys.argv[8])  # how many snapshots to save
cfl_const = float(sys.argv[9])

# OTHER OPTIONS

# true or false
nlinear=sys.argv[10]

fkernel=sys.argv[11]
#fkernel="zero_smooth_forcing"
#fkernel="zero_smooth_fourier_forcing"
#fkernel="zero_forcing"
#fkernel="white_forcing"
#fkernel="zero_white_forcing"

initial_value=sys.argv[12]
#initial_value="random_gaussian_smooth_initial"
#initial_value="zero_initial"
#initial_value="sine_initial"

scheme=sys.argv[13]
#scheme="euler"
#scheme="predcorr1"
#scheme="ETD"

saveformat=sys.argv[14]
# "fourier" or "real"




######## SET NUMERICAL PARAMETERS ###########
#############################################
#############################################
#############################################






path = os.path.abspath(os.path.join(os.path.abspath(''), os.pardir))
sys.path.append(path+'/scripts')

N = 2**BN
dx=1./N
L = Ltotal*Lrelative

# time parameters
#Ttotal = float(sys.argv[4]) # total time
dt = cfl_const*dx*dx  # time step size
NT = int(Ttotal/dt)  # number of time steps to reach stationary state

# calculated parameters
sqdt = m_sqrt(dt)
sqdx = m_sqrt(dx)
sqhdx = m_sqrt(.5*dx)
cte = 2.*pi/Ltotal
viscte = 4.*pi*pi*nu

# tolerance for zero values
#tol = 10.*np.finfo(float).eps
tol = 1e-14

#if c*Ttotal > .5*N:
#    raise Exception("Ttotal too big, larger than box size")

# final output will store NTsave snapshots of the whole spatial velocity field ...

skip = NT//NTsave
# ... and a few time evolutions (detailed in time, with Ntime instants, but only at a few positions)
Ntime  = 2**10
skit   = NT//Ntime
fewx   = [0,N//10,2*N//10,3*N//10,4*N//10,5*N//10,6*N//10,7*N//10,8*N//10,9*N//10]

# only works if workstation has the whole git folder
# get current commit short hash
#repo = Repo(search_parent_directories=True)
#git_hash = repo.head.object.hexsha[:10]

# print all parameters to aux file
all_params = {  "BN": BN, "Lrelative": Lrelative, "Ttotal": Ttotal,
                "N": N, "nu": nu, "Ltotal": Ltotal, "L": L, "sqeps": sqeps,
                "dx": dx, "NT": NT, "NTsave": NTsave,
                "viscte": viscte, "sqdx": sqdx,
                "cfl_const": cfl_const, "dt": dt,
                "nlinear": nlinear, "fkernel": fkernel,
                "initial_value": initial_value, "scheme": scheme,
                "saveformat": saveformat
            }

######## SAVE PARAMETERS WITH HASH ##########
#############################################
#############################################
#############################################

# to make sure all parameters are the same
# (with very high probability, since we do not consider the full hash)
# dumps means dump_String
dict_hash = sha224(dumps(all_params).encode('utf-8')).hexdigest()[:10]

# add these parameters to the list, they are ignored when
# calculating the hash
all_params["R"] = R
all_params["dict_hash"] = dict_hash

cwd = os.getcwd()
if R%1000 == 0: # print the parameter file
    with open(f'{cwd}/data/params_{dict_hash}.json','w',encoding="utf-8") as file:
        dump(all_params, file)





############ DEFINE FUNCTIONS ###############
#############################################
#############################################
#############################################






def getForce(f0):

    # new time, new forcing

    f0 = tc.normal(0.,1./sqdx,size=(N,),dtype=tc.cfloat)
    f0 = sqeps * kernel * fft(f0)

    return f0

""" In :    v1, field 1 in Fourier space
    Out:    v2=FFT( u1_x ), in Fourier space

    u1 = IFFT( v1 ), not needed explicitly
"""
def dudx(v1):

    v2 = 2.*pi*1j*K*v1

    return v2

""" In :    v1, field 1 in Fourier space
            v2, field 2 in Fourier space
    Out:    vf=FFT( u1 u2 ), in Fourier space
    Temp:
"""
def DealiasConvo(v1,v2):

    # Returns FFT(u1 u2) dealised where u1 = IFFT(v1) and u2 = IFFT(v2)

    M = N//2*3

    vpad1  = tc.zeros(M,dtype=tc.complex128)
    vpad2  = tc.zeros(M,dtype=tc.complex128)

    vpad1[0:N//2]   = v1[0:N//2]
    vpad1[M-N//2:M] = v1[N//2:]

    vpad2[0:N//2]   = v2[0:N//2]
    vpad2[M-N//2:M] = v2[N//2:]

    vpad1  = fft( fft(vpad1) * fft(vpad2) )
    vpad1 *= M/N

    vf = tc.zeros(N,dtype=tc.complex128)
    vf[0:N//2]   = vpad1[0:N//2]
    vf[N//2:]    = vpad1[M-N//2:M]

    return vf

""" In :    v0, velocity at previous instant
    Out:    v0, velocity at current instant
    Temp:   vp
            f0, force at current instant
"""
def BurgersEuler(v0,f0):

    vn  = tc.copy(v0)
    vn -= dt * viscte * K2 * v0 # viscous term
    vn += sqdt * f0 # force term
    if nlinear == 'True':
        vn -= dt * DealiasConvo( v0, dudx(v0) ) # nonlinear term

    return vn

""" In :    v0, velocity at previous instant
    Out:    vn, velocity at current instant
    Temp:   vp
            f0, force at current instant
"""
# Weak 1.0 Order Predictor-Corrector Method
# Kloeden-Platen PDFp.537, p. 502
def BurgersPredCorr1(v0,f0,vp):

    # predictor step
    vp  = tc.copy(v0)
    vp -= dt * viscte * K2 * v0 # viscous term
    vp += sqdt * f0 # force term
    if nlinear == 'True':
        vp -= dt * DealiasConvo( v0, dudx(v0) ) # nonlinear term

    vn  = tc.copy(v0)
    # half of corrector step
    vn -= .5 * dt * viscte * K2 * v0 # viscous term
    if nlinear == 'True':
        vn -= .5 * dt * DealiasConvo( v0, dudx(v0) ) # nonlinear term
    vn += sqdt * f0 # force term

    # other half of corrector step
    vn -= .5 * dt * viscte * K2 * vp # viscous term
    if nlinear == 'True':
        vn -= .5 * dt * DealiasConvo( vp, dudx(vp) ) # nonlinear term

    return vn

# Exponential Time Differencing (ETD) algorithm
# Forcing is treated together with nonlinear term
def BurgersETD(v0,f0):

    # nonlinear drift
    if nlinear == 'True':
        v0 -= dt * DealiasConvo( v0, dudx(v0) )
    v0 += sqdt * f0
    # linear drift, exponentiated
    v0 *= tc.exp(-viscte*dt*K2)

    # iron out imaginary part
    #v0 = fft( tc.real( fft(v0)) )

    return v0


# arrays
# simulation is usually run in Fourier space, even though end result is a real velocity field
# (hence Hermitian complex velocity fields)
vp = tc.zeros((N,),dtype=tc.complex128)
f0 = tc.zeros((N,),dtype=tc.complex128)

X = fftfreq(N) * Ltotal
K = fftfreq(N) * N
K2 = K*K

# choose initial condition
if   initial_value == "zero_initial":
    v0 = tc.zeros((N,),dtype=tc.complex128)
elif initial_value == "random_gaussian_smooth_initial":
    #v0 = fft( (1.-X*X/L**2)*exp(-.5*X**2/L/L) )
    #v0 = fft( exp(-.5*X**2/L**2) )
    raise NotImplementedError('Not done yet')
elif initial_value == "sine_initial":
    v0 = fft( tc.sin(2.*pi*X/Ltotal) )
else:
    raise ValueError('Invalid initial value string')


# CHOOSE forcing kernel
if   fkernel == "zero_smooth_forcing":
# gaussian, zero mode is zero, defined in real space
# make sure that zero-mode is exactly zero
# because of DFT, it has a small non-zero value
# apparently independent on N (a mistery to me)
    kernel = (1.-X*X/L**2)*tc.exp(-.5*X**2/L/L) # exponential correlation function
    kernel = tc.sqrt(fft(kernel))
    kernel[0] = 0.
elif fkernel == "smooth_forcing":
# make sure that zero-mode is exactly zero
# because of DFT, it has a small non-zero value
# apparently independent on N (a mistery to me)
    kernel = tc.exp(-.5*X**2/L/L) # exponential correlation function
    kernel = tc.sqrt(fft(kernel))
elif fkernel == "smooth_fourier_forcing":
    kernel = m_sqrt(2*pi*L*L)*tc.exp(-2*pi**2*L**2*K2) # exponential correlation function
    kernel = tc.sqrt(kernel)
elif fkernel == "zero_smooth_fourier_forcing":
# gaussian, zero mode is zero, defined in Fourier space
    kernel  = 4*m_sqrt(2*pi)*K2*L**3*pi**2
    kernel *= tc.exp(-2*pi**2*L**2*K2) # exponential correlation function
    kernel  = tc.sqrt(kernel)
elif fkernel == "zero_forcing":
# zero forcing
    kernel = tc.zeros((N,))
elif fkernel == "white_fourier_forcing":
# 3. white noise in real space, defined as flat in Fourier space
# all modes are populated, and normalized to have C_f(x=0) = 1. in real space
# to verify this, see tests/Test Fourier and Real Space Forcing Spectra.ipynb
    kernel = tc.ones(N)*m_sqrt(N)
    kernel = kernel / tc.sqrt(tc.count_nonzero(kernel))
elif fkernel == "zero_white_fourier_forcing":
# 4. white noise in real space, but zero mode is zero
# all modes are populated, and normalized to have C_f(x=0) = 1. in real space
# to verify this, see tests/Test Fourier and Real Space Forcing Spectra.ipynb
    kernel = tc.ones(N)*m_sqrt(N)
    kernel[0] = 0.
    kernel = kernel / tc.sqrt(tc.count_nonzero(kernel))
else:
    raise ValueError('Invalid forcing kernel string')

if saveformat != 'fourier' and saveformat != 'real':
    raise ValueError("Invalid save format string. Should be 'fourier' or 'real'")
if saveformat == 'real':
    raise NotImplementedError("Not implemented in the pytorch version")

# define output file
# and create HDF5 dataset
fname = f'{cwd}/data/burgers_R_{R:05d}_{dict_hash}'
f = File(fname, 'w')
dset = f.create_dataset(f'R_{R:05d}_{dict_hash}', (NT,N), dtype=complex128)

# initial steps of the fractional equation
# until it reaches a stationary state
for n in range(NTsave*skip):

    f0 = getForce(f0)

    # check that zero mode of force is zero
    assert tc.abs(f0[0]) < tol, f"zero mode of force is not zero, got: {f0[0]}"

    # check that forcing is purely real in real space
    assert tc.max(tc.imag(fft(f0))) < tol, f"forcing in real space is not real"

    if scheme=="euler":
        v0 = BurgersEuler(v0,f0)
    elif scheme=="predcorr1":
        v0 = BurgersPredCorr1(v0,f0,vp)
    elif scheme=="ETD":
        v0 = BurgersETD(v0,f0)

    # check that zero mode of velocity is zero
    assert tc.abs(v0[0]) < tol, f"zero mode of velocity is not zero, got: {v0[0]}"

    # check that velocity is purely real in real space
    assert tc.max(tc.imag(fft(v0))) < tol, f"velocity in real space is not real"

    #if saveformat == 'fourier':
    if (n+1)%skip == 0:
        dset[n//skip,:] = v0.detach().clone().numpy()

# close dataset
del dset
f.close()
