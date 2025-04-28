# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import scipy.constants
from skrf.media import Freespace
from skrf.frequency import Frequency

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

import pyWrapper as pW
#%%
nodalSourceLength = 0.01
VoltagePerMeterPeak = 1
resistanceLength = 0.01
resistancePerMeter = 1000
cableDiameter = 0.0001
# %%
fn = "NodalSourceTest.fdtd.json"
solver = pW.FDTD(input_filename=fn, path_to_exe=SEMBA_EXE)
# %%
sourcePointProbe = pW.Probe(solver.getSolvedProbeFilenames("Point probe bottom")[0])
sourceBulkProbe = pW.Probe(solver.getSolvedProbeFilenames("Bulk probe Nodal Source")[0])
resistancePointProbe = pW.Probe(solver.getSolvedProbeFilenames("Point probe top")[0])
resistanceBulkProbe = pW.Probe(solver.getSolvedProbeFilenames("Bulk probe Resistance")[0])
excitation = pW.ExcitationFile(excitation_filename=solver.getExcitationFile("predefinedExcitation")[0])
#%%
sourcePointProbe.data['voltage'] = sourcePointProbe.data['field'] * nodalSourceLength

print("Max field: ", sourcePointProbe.data['field'].max())
print("Max voltage: ", sourcePointProbe.data['voltage'].max())
print("Max Current Readed: ", sourceBulkProbe.data['current'].max())
#%%
resistancePointProbe.data['voltage'] = resistancePointProbe.data['field'] * resistanceLength
resistancePointProbe.data['current'] = resistancePointProbe.data['voltage'] / (resistancePerMeter/resistanceLength)
print("Max field: ", resistancePointProbe.data['field'].max())
print("Max voltage: ", resistancePointProbe.data['voltage'].max())
print("Max current: ", resistancePointProbe.data['current'].max())
plt.plot(resistancePointProbe.data['current'], label='Expected Current')
plt.plot(resistanceBulkProbe.data['current'], label='Readed Current')
plt.legend()
plt.show()
print(resistanceBulkProbe.data['current'].max())
print(resistancePointProbe.data['current'].max())

#assert np.allclose(resistancePointProbe.data['current'], resistanceBulkProbe.data['current'])
#%%













resistancePointProbe.data['voltage'] = resistancePointProbe.data['field'] * resistanceLength
resistancePointProbe.data['current'] = resistancePointProbe.data['voltage'] / (resistancePerMeter/resistanceLength)
# %%
fieldResistanceMax = resistancePointProbe.data['field'].abs().max()
timeFieldResistanceMax = resistancePointProbe.data.loc[resistancePointProbe.data['field'] == fieldResistanceMax, 'time'].iloc[0]
voltageFieldResistanceMax = resistancePointProbe.data.loc[resistancePointProbe.data['time'] == timeFieldResistanceMax, 'voltage'].iloc[0]
currentFieldResistanceMax = resistancePointProbe.data.loc[resistancePointProbe.data['time'] == timeFieldResistanceMax, 'current'].iloc[0]
print("Resistance Values:")
print("Max field: ", fieldResistanceMax)
print("Time at max field: ", timeFieldResistanceMax)
print("Voltage at max field: ", voltageFieldResistanceMax)
print("Current at max field: ", currentFieldResistanceMax)
#%%
sourcePointProbe.data['voltage'] = sourcePointProbe.data['field'] * nodalSourceLength

fieldSourceMax = sourcePointProbe.data['field'].abs().max()
timeFieldSourceMax = sourcePointProbe.data.loc[sourcePointProbe.data['field'] == fieldSourceMax, 'time'].iloc[0]
voltageFieldSourceMax = sourcePointProbe.data.loc[sourcePointProbe.data['time'] == timeFieldSourceMax, 'voltage'].iloc[0]

print("Source Values:")
print("Max field: ", fieldSourceMax)
print("Time at max field: ", timeFieldSourceMax)
print("Voltage at max field: ", voltageFieldSourceMax)

#%%
voltageDelay = timeFieldResistanceMax - timeFieldSourceMax
lightDelayTime = 0.120 / 3e08

print(voltageDelay)
print(lightDelayTime)

#%%
resistancePointProbe.data['field_normalized'] = resistancePointProbe.data['field'] / resistancePointProbe.data['field'].abs().max()
plt.plot(resistancePointProbe.data['time'], resistancePointProbe.data['field_normalized'], label="Probe")
plt.plot(excitation.data['time'], excitation.data['value'], label = "Excitation")
plt.legend()
plt.show()

# %%
dt = resistancePointProbe.data['time'].diff().fillna(0)
dq = resistancePointProbe.data['current'] * dt
resistancePointProbe.data['charge'] = dq.cumsum()
total_charge = resistancePointProbe.data.loc[resistancePointProbe.data['time'] == timeFieldSourceMax, 'time'].iloc[0]
print("Total accumulated charge on the resistance:", total_charge)
print("Capacity: ", total_charge/voltageFieldSourceMax)


# %%
#plt.plot(resistancePointProbe['voltage'], label='voltage')
plt.plot(resistancePointProbe['current'], label='current')
plt.legend()
plt.show()
# %%
plt.plot(resistanceBulkProbe.data['current'])
plt.show()
# %%
plt.plot(resistanceBulkProbe.data['current'], label="Resisitance")
plt.plot(-1*sourceBulkProbe.data['current'], label = "Nodal")
plt.legend()
plt.show()
# %%
plt.plot(resistanceBulkProbe.data['current'], label="Resisitance. Max: {}".format(resistanceBulkProbe.data['current'].max()))
plt.plot(-1*sourceBulkProbe.data['current'], label = "Nodal. Max: {}".format(-1*sourceBulkProbe.data['current'].max()))
plt.legend()
plt.xlim(4000, 10000)
plt.show()
# %%
