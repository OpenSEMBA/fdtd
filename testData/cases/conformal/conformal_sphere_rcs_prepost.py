# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import scipy.constants
from scipy.special import hankel2 as h
from scipy.special import h2vp as hp

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *


# %% Generate excitation and visualize
dt = 1e-12
w0 = 0.1e-9    # ~ 2 GHz bandwidth
# t0 = 10*w0
t0 = 1e-9
t = np.arange(0, t0+20*w0, dt)
f = np.exp( -np.power(t-t0,2)/ w0**2 )

plt.figure()
plt.plot(t,f)

data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = f
np.savetxt('gauss.exc', data)

# %% Prepare solver
solver = FDTD(input_filename = 'conformal_sphere_1mm_rcs_delta.fdtd.json', path_to_exe=SEMBA_EXE)

solver.cleanUp()
solver.run()

# %% analytical result

def h_sph(n,z):
    return np.sqrt(0.5*np.pi*z)*h(n+0.5,z)

def h_sph_der(n,z):
    return 0.25*np.pi*(0.5*np.pi*z)**(-0.5)*h(n+0.5,z) + np.sqrt(0.5*np.pi*z)*hp(n+0.5,z)

def ha(n,z):
    return (1j)**(n+1)*np.exp(-1j*z)
def had(n,z):
    return (1j)**n*np.exp(-1j*z)

a = 0.01 #radius
f = np.linspace(3e8,2.4e11,1000)
l = 3.0e8/f
beta = 2*np.pi*f/3.0e8

rcs = 0.0
for n in range(1,50):
    rcs = rcs + (-1)**n*(2*n+1)/(h_sph_der(n,beta*a)*h_sph(n,beta*a))
    # rcs = rcs + ((-1)**n*(2*n+1)/(ha(n,beta*a)*had(n, beta*a)))**2
rcs = np.abs(rcs)**2*(l**2/(4*np.pi))


# %% plot results 

far = Probe(solver.getSolvedProbeFilenames("Far")[0])
ra0_0 = far.data.loc[(far.data['Theta'] == 0) & (far.data['Phi'] == 0), 'rcs_arit']
rg0_0 = far.data.loc[(far.data['Theta'] == 0) & (far.data['Phi'] == 0), 'rcs_geom']
f0_0  = far.data.loc[(far.data['Theta'] == 0) & (far.data['Phi'] == 0), 'freq']

fieldx = Probe(solver.getSolvedProbeFilenames("electric_field_point_probe_Ex")[0])
fieldy = Probe(solver.getSolvedProbeFilenames("electric_field_point_probe_Ey")[0])
fieldz = Probe(solver.getSolvedProbeFilenames("electric_field_point_probe_Ez")[0])

plt.plot(fieldx['time'], fieldx['incident'], label='incident x')
plt.plot(fieldx['time'], fieldy['incident'], label='incident y')
plt.plot(fieldx['time'], fieldz['incident'], label='incident z')
plt.legend()
plt.show()
plt.plot(fieldx['time'], fieldx['field'], label='field x')
plt.plot(fieldx['time'], fieldy['field'], label='field y')
plt.plot(fieldx['time'], fieldz['field'], label='field z')
plt.legend()
plt.show()
# freqs = far.data['freq'].unique()
# rcs_sum = np.zeros(0)
# for freq in freqs:
#     rcs_sum = np.append(rcs_sum, far.data.loc[(far.data['freq'] == freq) & (far.data['Theta'] == 10.0) ].rcs_arit.sum())

# #normalize by number of phi values?
# plt.plot(far['freq'], far['rcs_arit'], '-', label='simulation')
# plt.plot(f0_0, rg0_0/(np.pi*(1000)**2), '.', label='rcs arit')
plt.plot(f0_0, ra0_0/(np.pi*(1000)**2), '.', label='rcs geom')
plt.plot(f, rcs/(np.pi*a**2), '-.', label='analytical')
# plt.plot(a/l, rcs/(np.pi*a**2), '-.', label='analytical')
plt.xlabel('f [Hz]')
plt.ylabel(r'rcs/$\pi*a^2$')
plt.title('Sphere radar cross section')
# plt.hlines(1,0.0,5e8, colors ='red', linestyles='--')
plt.legend()


# %% Postprocess
