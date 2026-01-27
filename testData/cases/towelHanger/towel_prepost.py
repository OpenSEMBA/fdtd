# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt


import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *

def createWire(id, r, rpul = 0.0, lpul=0.0):
    return {"id":id,
            "type": "wire",
            "radius": r, 
            "resistancePerMeter": rpul, 
            "inductancePerMeter": lpul
            }

def createUnshieldedWire(id, lpul, cpul):
    return {         
        "id": id,
        "type": "unshieldedMultiwire",
        "inductancePerMeter" :  [[lpul]],
        "capacitancePerMeter" : [[cpul]],
        "resistancePerMeter": [0.0],
        "conductancePerMeter": [0.0]
        }

#####################################################
# %% Run solver
fn = 'towelHanger.fdtd.json'

# exe = os.path.join(os.getcwd(), 'build-dbg-nomtln', 'bin', 'semba-fdtd')
solver = FDTD(input_filename = fn, path_to_exe='../../../build-dbg-nomtln/bin/semba-fdtd')
solver['materials'][0] = createWire(id = 1, r = 0.1e-3)
solver.cleanUp()
solver.run()

# %% Run solver
p_wire_0 = Probe(solver.getSolvedProbeFilenames("wire_start")[0])
p_wire_1 = Probe(solver.getSolvedProbeFilenames("wire_mid")[0])
p_wire_2 = Probe(solver.getSolvedProbeFilenames("wire_end")[0])

# %% Run solver

solver = FDTD(input_filename = fn, path_to_exe='../../../build-dbg/bin/semba-fdtd')
solver['materials'][0] = createUnshieldedWire(id = 1, lpul = 6.5183032590978384e-07, cpul = 1.7046017451862063e-11)        
solver.cleanUp()
solver.run()

p_mwire_0 = Probe(solver.getSolvedProbeFilenames("wire_start")[0])
p_mwire_1 = Probe(solver.getSolvedProbeFilenames("wire_mid")[0])
p_mwire_2 = Probe(solver.getSolvedProbeFilenames("wire_end")[0])
#####################################################
# %% Postprocess
OUTPUTS_FOLDER = '../../../testData/outputs/'
expected_0 = Probe(OUTPUTS_FOLDER+'towelHanger.fdtd_wire_start_Wz_27_25_30_s1.dat')
expected_1 = Probe(OUTPUTS_FOLDER+'towelHanger.fdtd_wire_mid_Wx_35_25_32_s5.dat')
expected_2 = Probe(OUTPUTS_FOLDER+'towelHanger.fdtd_wire_end_Wz_43_25_30_s4.dat')


plt.figure()
plt.plot(expected_0['time']*1e9, expected_0['current'],'-', label='expected')
plt.plot(p_wire_0['time']*1e9, p_wire_0['current'],'-', label='wire')
plt.plot(p_mwire_0['time']*1e9, -p_mwire_0['current_0'],'-', label='multiwire')
plt.xlabel('Time (ns)')
plt.ylabel('Current (A)')
plt.grid(which='both')
plt.legend()

plt.figure()
plt.plot(expected_1['time']*1e9, expected_1['current'],'-', label='expected')
plt.plot(p_wire_1['time']*1e9, p_wire_1['current'],'-', label='wire')
plt.plot(p_mwire_1['time']*1e9, -p_mwire_1['current_0'],'-', label='multiwire')
plt.xlabel('Time (ns)')
plt.ylabel('Current (A)')
plt.grid(which='both')
plt.legend()

plt.figure()
plt.plot(expected_2['time']*1e9, expected_2['current'],'-', label='expected')
plt.plot(p_wire_2['time']*1e9, p_wire_2['current'],'-', label='wire')
plt.plot(p_mwire_2['time']*1e9, p_mwire_2['current_0'],'-', label='multiwire')
plt.xlabel('Time (ns)')
plt.ylabel('Current (A)')
plt.grid(which='both')
plt.legend()


# %%
print (np.allclose(-expected_0['current'], p_mwire_0['current_0'],  rtol=5e-3, atol=5e-4))
print (np.allclose(-expected_1['current'], p_mwire_1['current_0'],  rtol=5e-3, atol=5e-4))
print (np.allclose(-expected_2['current'], -p_mwire_2['current_0'],  rtol=5e-3, atol=5e-4))
# %%

mu0 = 4*np.pi*1e-7
kl = np.linspace(0,0.1, 1000)
plt.figure()
plt.loglog(kl, mu0/(mu0 + 100*mu0*kl), label= 'r=100')
plt.loglog(kl, mu0/(mu0 + 1000*mu0*kl), label= 'r=1000')
plt.loglog(kl, mu0/(mu0 + 10000*mu0*kl), label= 'r=10000')
plt.legend()
plt.plot()

# %%
