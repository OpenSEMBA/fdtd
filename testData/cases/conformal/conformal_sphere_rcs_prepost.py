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
solver = FDTD(input_filename = 'conformal_str_sphere_rcs.fdtd.json', path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

far_str = Probe(solver.getSolvedProbeFilenames("n2f")[0])
ra90_0_str = far_str.data.loc[(far_str.data['Theta'] == 90.0) & (far_str.data['Phi'] ==  0.0), 'rcs_arit']
rg90_0_str = far_str.data.loc[(far_str.data['Theta'] == 90.0) & (far_str.data['Phi'] ==  0.0), 'rcs_geom']
f90_0_str  = far_str.data.loc[(far_str.data['Theta'] == 90.0) & (far_str.data['Phi'] ==  0.0), 'freq']



# %% Slanted
solver = FDTD(input_filename = 'slanted_sphere_rcs.fdtd.json', path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

far_sl = Probe(solver.getSolvedProbeFilenames("n2f")[0])
ra90_0_sl = far_sl.data.loc[(far_sl.data['Theta'] == 90.0) & (far_sl.data['Phi'] ==  0.0), 'rcs_arit']
rg90_0_sl = far_sl.data.loc[(far_sl.data['Theta'] == 90.0) & (far_sl.data['Phi'] ==  0.0), 'rcs_geom']
f90_0_sl  = far_sl.data.loc[(far_sl.data['Theta'] == 90.0) & (far_sl.data['Phi'] ==  0.0), 'freq']


# %% Conformal
solver = FDTD(input_filename = 'conformal_fL_0.15_sphere_rcs.fdtd.json', path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

far_conf = Probe(solver.getSolvedProbeFilenames("n2f")[0])
ra90_0_conf = far_conf.data.loc[(far_conf.data['Theta'] == 90.0) & (far_conf.data['Phi'] ==  0.0), 'rcs_arit']
rg90_0_conf = far_conf.data.loc[(far_conf.data['Theta'] == 90.0) & (far_conf.data['Phi'] ==  0.0), 'rcs_geom']
f90_0_conf  = far_conf.data.loc[(far_conf.data['Theta'] == 90.0) & (far_conf.data['Phi'] ==  0.0), 'freq']


# %% Conformal finer
solver = FDTD(input_filename = 'conformal_fL_sphere_rcs.fdtd.json', path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

far_conf_fine = Probe(solver.getSolvedProbeFilenames("n2f")[0])
ra90_0_conf_fine = far_conf_fine.data.loc[(far_conf_fine.data['Theta'] == 90.0) & (far_conf_fine.data['Phi'] ==  0.0), 'rcs_arit']
rg90_0_conf_fine = far_conf_fine.data.loc[(far_conf_fine.data['Theta'] == 90.0) & (far_conf_fine.data['Phi'] ==  0.0), 'rcs_geom']
f90_0_conf_fine  = far_conf_fine.data.loc[(far_conf_fine.data['Theta'] == 90.0) & (far_conf_fine.data['Phi'] ==  0.0), 'freq']

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
plt.loglog(f90_0_str, rg90_0_str/(np.pi*(a)**2), '--', label='Structured')
plt.loglog(f90_0_sl, rg90_0_sl/(np.pi*(a)**2), '--', label='Slanted')
plt.loglog(f90_0_conf, rg90_0_conf/(np.pi*(a)**2), '.', label='Conformal fL = 0.15; eP = 5, ')
plt.loglog(f90_0_conf, rg90_0_conf/(np.pi*(a)**2), '.', label='Conformal fL = 0.005; eP = 50')
plt.xlabel('f [Hz]')
plt.ylabel(r'rcs/$\pi*a^2$')
plt.title('Sphere radar cross section')
plt.hlines(1,0.0,7e8, colors ='red', linestyles='--')
plt.legend()

# %%
rcs_interp = np.interp(
    f,
    f90_0_conf,
    rg90_0_conf
)

assert np.allclose(rcs[5:150], rcs_interp[5:150], rtol=0.25, atol=0.15)
