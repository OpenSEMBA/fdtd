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


# %% Structured
solver = FDTD(input_filename = 'conformal_surface_sphere.fdtd.json', path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

# %% Structured
far = Probe(solver.getSolvedProbeFilenames("n2f")[0])
ra90_0 = far.data.loc[(far.data['Theta'] == 90.0) & (far.data['Phi'] ==  0.0), 'rcs_arit']
rg90_0 = far.data.loc[(far.data['Theta'] == 90.0) & (far.data['Phi'] ==  0.0), 'rcs_geom']
f90_0  = far.data.loc[(far.data['Theta'] == 90.0) & (far.data['Phi'] ==  0.0), 'freq']

ex = Probe(solver.getSolvedProbeFilenames("indside_Ex")[0])
ey = Probe(solver.getSolvedProbeFilenames("indside_Ey")[0])
ez = Probe(solver.getSolvedProbeFilenames("indside_Ez")[0])

# %% analytical result

def h_sph(n,z):
    return np.sqrt(0.5*np.pi*z)*h(n+0.5,z)

def h_sph_der(n,z):
    return 0.25*np.pi*(0.5*np.pi*z)**(-0.5)*h(n+0.5,z) + np.sqrt(0.5*np.pi*z)*hp(n+0.5,z)

def ha(n,z):
    return (1j)**(n+1)*np.exp(-1j*z)
def had(n,z):
    return (1j)**n*np.exp(-1j*z)

a = 0.5 #radius
f = np.linspace(1e7,6e8,201)
l = 3.0e8/f
beta = 2*np.pi*f/3.0e8

rcs = 0.0
for n in range(1,50):
    rcs = rcs + (-1)**n*(2*n+1)/(h_sph_der(n,beta*a)*h_sph(n,beta*a))
rcs = np.abs(rcs)**2*(l**2/(4*np.pi))


# %% plot results 
plt.loglog(f, rcs/(np.pi*a**2), '.', label='analytical')
plt.loglog(f90_0, rg90_0/(np.pi*(a)**2), '--', label='Structured')
plt.xlabel('f [Hz]')
plt.ylabel(r'rcs/$\pi*a^2$')
plt.title('Sphere radar cross section')
plt.hlines(1,0.0,7e8, colors ='red', linestyles='--')
plt.legend()

# %% plot results 
plt.plot(ex['time']*1e6, ex['field'], label=r'$E_x$')
plt.plot(ey['time']*1e6, ey['field'], label=r'$E_y$')
plt.plot(ez['time']*1e6, ez['field'], label=r'$E_z$')
plt.xlabel(r't [$\mu$s]')
plt.ylabel(r'|E|')
plt.title('Field inside hollow PEC sphere')
plt.legend()

