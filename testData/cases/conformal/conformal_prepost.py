# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import scipy.constants

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
solver = FDTD(input_filename = 'conformal.fdtd.json', path_to_exe=SEMBA_EXE)

dx = solver['mesh']['grid']['steps']['x'][0]

solver['materialAssociations'][0]['elementIds'] = [4]
solver['mesh']['elements'][3]['intervals'] = [[[0,0,4],[2,2,4]]]
solver.cleanUp()

solver.run()
front = Probe(solver.getSolvedProbeFilenames("front")[0])
t = front['time']
t4 = t[front['field'].argmin()]


solver['materialAssociations'][0]['elementIds'] = [4]
solver['mesh']['elements'][3]['intervals'] = [[[0,0,5],[2,2,5]]]
solver.cleanUp()

solver.run()
front = Probe(solver.getSolvedProbeFilenames("front")[0])
t = front['time']
t5 = t[front['field'].argmin()]

z = np.zeros(11)
td = np.zeros(11)
# for (i,delta) in enumerate([0.1]):
solver['materialAssociations'][0]['elementIds'] = [5]
for i in range(11):
    solver['mesh']['coordinates'][6]["relativePosition"][2]   = 4.0 + i*0.1
    solver['mesh']['coordinates'][7]["relativePosition"][2]   = 4.0 + i*0.1
    solver['mesh']['coordinates'][8]["relativePosition"][2]   = 4.0 + i*0.1
    solver['mesh']['coordinates'][9]["relativePosition"][2]  = 4.0  + i*0.1
    solver['mesh']['coordinates'][15]["relativePosition"][2]  = 4.0 + i*0.1
    solver['mesh']['coordinates'][16]["relativePosition"][2]  = 4.0 + i*0.1
    solver['mesh']['coordinates'][17]["relativePosition"][2]  = 4.0 + i*0.1
    solver['mesh']['coordinates'][18]["relativePosition"][2]  = 4.0 + i*0.1
    solver['mesh']['coordinates'][19]["relativePosition"][2]  = 4.0 + i*0.1

    solver.cleanUp()
    solver.run()
    assert solver.hasFinishedSuccessfully()
    front = Probe(solver.getSolvedProbeFilenames("front")[0])
    t = front['time']
    inc = front['incident']
    back = front['field']
    z[i] = i*0.1
    td[i] = t[back.argmin()]


# %% Postprocess
front = Probe(solver.getSolvedProbeFilenames("front")[0])
back = Probe(solver.getSolvedProbeFilenames("back")[0])

plt.plot(front['time'], front['incident'], 'k-', label='front incident') 
plt.plot(front['time'], front['field'], 'b.', label='front field') 
plt.legend()
plt.xlim(0,6e-9)


# %%
plt.plot(z,(td-t4)*1e9, '*r-', label=r'conformal cell 4 < z < 5') 
plt.axhline(y = 0, color = 'b', label = 'PEC surface at z = 4')
plt.axhline(y = (t5-t4)*1e9, color = 'g', label = 'PEC surface at z = 5')
plt.xlabel(r'conformal length $\Delta$ [cells]')
plt.ylabel(r'$\Delta$ t [ns]')
plt.title('Delay in reflection wrt to a PEC box in z = 4')
plt.legend()

# %%
for i in range(10):
    print(td[i], ' ', t4 + 2*(i*1.0/10)*dx/3e8)
# %%
