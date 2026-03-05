# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
from os import environ as env
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build-rls/bin/semba-fdtd'

from pyWrapper import *
env["SPICE_SCRIPTS"] = "./"

# %%


# %% Generate excitation and visualize
dt = 1.0e-12
t = np.arange(0, 3e-7, dt)
A = 45
freq = 80e6
f = A*np.sin(2*np.pi*freq*t)
plt.figure()
plt.plot(t,f)

data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = f
np.savetxt('sin_80M_45Vm.exc', data) 

# %% Generate excitation and visualize
dt = 0.1e-13
w0 = 0.05e-8 
t0 = 20*w0
t = np.arange(0, t0+20*w0, dt)
A = 300
f = A*np.exp( -np.power(t-t0,2)/ w0**2 )

plt.figure()
plt.plot(t,f)

data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = f
np.savetxt('gauss_300Vm_80M_500M.exc', data) 

# fq = fftfreq(len(t))/dt
# F = fft(f)
# F = F[fq>=0]
# fq = fq[fq>=0]
# Fmax = np.max(np.abs(F))
# Fcut = Fmax / np.sqrt(2.0)
# idx = (np.abs(np.abs(F) - Fcut)).argmin()

# # print(Fmax, Fcut)

# plt.figure()
# plt.plot(fq, np.abs(F),'.-')
# plt.axvline(x = fq[idx], color = 'r')
# plt.xlim(80e6,1000e6)
# plt.ylim(1e0, 1e1)
# # plt.yscale('log')

# plt.grid()

# %% Prepare solver case debug
solver = FDTD(input_filename = 'wireCoupling_debug.fdtd.json', path_to_exe=SEMBA_EXE)
# solver["materialAssociations"][1]["initialTerminalId"] = 4
# solver["materials"][5]["relativePermittivity"] = 2.5
# solver["materials"][6]["relativePermittivity"] = 1.3
solver.cleanUp()
solver.run() 

v_input =  Probe(solver.getSolvedProbeFilenames("circuit_probe_wire")[0])
v_output =  Probe(solver.getSolvedProbeFilenames("circuit_probe_segment")[0])


# %% Prepare solver case A
solver = FDTD(input_filename = 'wireCoupling.fdtd.json', path_to_exe=SEMBA_EXE)
# solver["mesh"]["elements"][0]["coordinateIds"] = [6,1,4,2,3]
# solver["mesh"]["elements"][10]["coordinateIds"] = [5,6]
# solver["mesh"]["elements"][11]["coordinateIds"] = [6]
# solver["materialAssociations"][1]["initialTerminalId"] = 4
# solver["materials"][5]["relativePermittivity"] = 2.5
# solver["materials"][6]["relativePermittivity"] = 1.3
solver.cleanUp()
solver.run() 

v_input =  Probe(solver.getSolvedProbeFilenames("circuit_probe_wire")[0])
v_output =  Probe(solver.getSolvedProbeFilenames("circuit_probe_segment")[0])

#%% sava data
data = np.zeros((len(v_input['time']), 3))
data[:,0] = v_input['time']
data[:,1] = v_input['voltage_0']
data[:,2] = v_output['voltage_0']
np.savetxt('../../outputs/wireCoupling/caseA_exc_sin_80M_45Vm.dat', data)

# %% Plot time domain
offset = 0.0375
plt.figure()
plt.plot(v_input['time'], v_input['voltage_0'],'-', label='input, no plane')
plt.plot(v_output['time'], v_output['voltage_0']-offset,'--', label='output, no plane')
plt.grid(which='both')
plt.xlabel('time [s]')
plt.ylabel('V [V]')
# plt.xlim(0,0.8e-7)
# plt.ylim(-3,3)
plt.legend()

#%% 
# solver = FDTD(input_filename = 'wireCouplingH.fdtd.json', path_to_exe=SEMBA_EXE)
# solver["mesh"]["elements"][0]["coordinateIds"] = [6,1,4,2,3]
# solver["mesh"]["elements"][10]["coordinateIds"] = [5,6]
# solver["mesh"]["elements"][11]["coordinateIds"] = [6]
# solver["materialAssociations"][1]["initialTerminalId"] = 2
# solver["materials"][5]["relativePermittivity"] = 2.5
# solver["materials"][6]["relativePermittivity"] = 1.3
# solver.cleanUp()
# solver.run() 
# v_inputH =  Probe(solver.getSolvedProbeFilenames("circuit_probe_wire")[0])
# v_outputH =  Probe(solver.getSolvedProbeFilenames("circuit_probe_segment")[0])

# solver = FDTD(input_filename = 'wireCouplingV.fdtd.json', path_to_exe=SEMBA_EXE)
# solver["mesh"]["elements"][0]["coordinateIds"] = [6,1,4,2,3]
# solver["mesh"]["elements"][10]["coordinateIds"] = [5,6]
# solver["mesh"]["elements"][11]["coordinateIds"] = [6]
# solver["materialAssociations"][1]["initialTerminalId"] = 2
# solver["materials"][5]["relativePermittivity"] = 2.5
# solver["materials"][6]["relativePermittivity"] = 1.3
# solver.cleanUp()
# solver.run() 
# v_inputV =  Probe(solver.getSolvedProbeFilenames("circuit_probe_wire")[0])
# v_outputV =  Probe(solver.getSolvedProbeFilenames("circuit_probe_segment")[0])

# %% Prepare solver case B
solver = FDTD(input_filename = 'wireCoupling.fdtd.json', path_to_exe=SEMBA_EXE)
solver["mesh"]["elements"][0]["coordinateIds"] = [1,4,2,3]
solver["mesh"]["elements"][10]["coordinateIds"] = [5,1]
solver["mesh"]["elements"][11]["coordinateIds"] = [1]
solver["materialAssociations"][1]["initialTerminalId"] = 2
solver.cleanUp()
solver.run() 
v_input =  Probe(solver.getSolvedProbeFilenames("circuit_probe_wire")[0])
v_output =  Probe(solver.getSolvedProbeFilenames("circuit_probe_segment")[0])

solver = FDTD(input_filename = 'wireCouplingH.fdtd.json', path_to_exe=SEMBA_EXE)
solver["mesh"]["elements"][0]["coordinateIds"] = [1,4,2,3]
solver["mesh"]["elements"][10]["coordinateIds"] = [5,1]
solver["mesh"]["elements"][11]["coordinateIds"] = [1]
solver["materialAssociations"][1]["initialTerminalId"] = 2
solver.cleanUp()
solver.run() 
v_inputH =  Probe(solver.getSolvedProbeFilenames("circuit_probe_wire")[0])
v_outputH =  Probe(solver.getSolvedProbeFilenames("circuit_probe_segment")[0])

solver = FDTD(input_filename = 'wireCouplingV.fdtd.json', path_to_exe=SEMBA_EXE)
solver["mesh"]["elements"][0]["coordinateIds"] = [1,4,2,3]
solver["mesh"]["elements"][10]["coordinateIds"] = [5,1]
solver["mesh"]["elements"][11]["coordinateIds"] = [1]
solver["materialAssociations"][1]["initialTerminalId"] = 2
solver.cleanUp()
solver.run() 
v_inputV =  Probe(solver.getSolvedProbeFilenames("circuit_probe_wire")[0])
v_outputV =  Probe(solver.getSolvedProbeFilenames("circuit_probe_segment")[0])

# %% Prepare solver case C
solver = FDTD(input_filename = 'wireCoupling.fdtd.json', path_to_exe=SEMBA_EXE)
solver["mesh"]["elements"][0]["coordinateIds"] = [6,1,4,2,3]
solver["mesh"]["elements"][10]["coordinateIds"] = [5,6]
solver["mesh"]["elements"][11]["coordinateIds"] = [6]
solver["materialAssociations"][1]["initialTerminalId"] = 3
solver.cleanUp()
solver.run() 
v_input =  Probe(solver.getSolvedProbeFilenames("circuit_probe_wire")[0])
v_output =  Probe(solver.getSolvedProbeFilenames("circuit_probe_segment")[0])

solver = FDTD(input_filename = 'wireCouplingH.fdtd.json', path_to_exe=SEMBA_EXE)
solver["mesh"]["elements"][0]["coordinateIds"] = [6,1,4,2,3]
solver["mesh"]["elements"][10]["coordinateIds"] = [5,6]
solver["mesh"]["elements"][11]["coordinateIds"] = [6]
solver["materialAssociations"][1]["initialTerminalId"] = 3
solver.cleanUp()
solver.run() 
v_inputH =  Probe(solver.getSolvedProbeFilenames("circuit_probe_wire")[0])
v_outputH =  Probe(solver.getSolvedProbeFilenames("circuit_probe_segment")[0])

solver = FDTD(input_filename = 'wireCouplingV.fdtd.json', path_to_exe=SEMBA_EXE)
solver["mesh"]["elements"][0]["coordinateIds"] = [6,1,4,2,3]
solver["mesh"]["elements"][10]["coordinateIds"] = [5,6]
solver["mesh"]["elements"][11]["coordinateIds"] = [6]
solver["materialAssociations"][1]["initialTerminalId"] = 3
solver.cleanUp()
solver.run() 
v_inputV =  Probe(solver.getSolvedProbeFilenames("circuit_probe_wire")[0])
v_outputV =  Probe(solver.getSolvedProbeFilenames("circuit_probe_segment")[0])
# %% Prepare solver case D
solver = FDTD(input_filename = 'wireCoupling.fdtd.json', path_to_exe=SEMBA_EXE)
solver["mesh"]["elements"][0]["coordinateIds"] = [1,4,2,3]
solver["mesh"]["elements"][10]["coordinateIds"] = [5,1]
solver["mesh"]["elements"][11]["coordinateIds"] = [1]
solver["materialAssociations"][1]["initialTerminalId"] = 3
solver.cleanUp()
solver.run() 
v_input =  Probe(solver.getSolvedProbeFilenames("circuit_probe_wire")[0])
v_output =  Probe(solver.getSolvedProbeFilenames("circuit_probe_segment")[0])

solver = FDTD(input_filename = 'wireCouplingH.fdtd.json', path_to_exe=SEMBA_EXE)
solver["mesh"]["elements"][0]["coordinateIds"] = [1,4,2,3]
solver["mesh"]["elements"][10]["coordinateIds"] = [5,1]
solver["mesh"]["elements"][11]["coordinateIds"] = [1]
solver["materialAssociations"][1]["initialTerminalId"] = 3
solver.cleanUp()
solver.run() 
v_inputH =  Probe(solver.getSolvedProbeFilenames("circuit_probe_wire")[0])
v_outputH =  Probe(solver.getSolvedProbeFilenames("circuit_probe_segment")[0])

solver = FDTD(input_filename = 'wireCouplingV.fdtd.json', path_to_exe=SEMBA_EXE)
solver["mesh"]["elements"][0]["coordinateIds"] = [1,4,2,3]
solver["mesh"]["elements"][10]["coordinateIds"] = [5,1]
solver["mesh"]["elements"][11]["coordinateIds"] = [1]
solver["materialAssociations"][1]["initialTerminalId"] = 3
solver.cleanUp()
solver.run() 
v_inputV =  Probe(solver.getSolvedProbeFilenames("circuit_probe_wire")[0])
v_outputV =  Probe(solver.getSolvedProbeFilenames("circuit_probe_segment")[0])



#%%
plt.figure()
plt.plot(v_inputH['time'], v_inputH['voltage_0'],'.-', label='input, H plane')
plt.plot(v_outputH['time'], v_outputH['voltage_0'],'--', label='output, H plane')
plt.grid(which='both')
plt.xlabel('time [s]')
plt.ylabel('V [V]')
plt.legend()

plt.figure()
plt.plot(v_inputV['time'], v_inputV['voltage_0'],'.-', label='input, V plane')
plt.plot(v_outputV['time'], v_outputV['voltage_0'],'--', label='output, V plane')
plt.grid(which='both')
plt.xlabel('time [s]')
plt.ylabel('V [V]')
plt.legend()

plt.show()
# # %% Prepare solver no plane

# solver["materialAssociations"][1]["initialTerminalId"] = 4
# solver.cleanUp()
# solver.run() 
# v_input_50 =  Probe(solver.getSolvedProbeFilenames("circuit_probe_wire")[0])
# v_output_50 =  Probe(solver.getSolvedProbeFilenames("circuit_probe_segment")[0])

# # %% output, no plane
# v_output_1M_r =  Probe(solver.getSolvedProbeFilenames("r_probe")[0])
# v_output_50_r =  Probe(solver.getSolvedProbeFilenames("r_probe")[0])

# plt.figure()
# plt.plot(v_output_1M['time'], v_output_1M['voltage_0'],'.-', label='output 1M, no plane')
# plt.plot(v_output_50['time'], v_output_50['voltage_0'],'--', label='output 50, no plane')
# plt.plot(v_output_1M_r['time'], v_output_1M_r['voltage_0'],'.-', label='output 1M r, no plane')
# plt.plot(v_output_50_r['time'], v_output_50_r['voltage_0'],'--', label='output 50 r, no plane')
# plt.grid(which='both')
# plt.xlabel('time [s]')
# plt.ylabel('V [V]')
# # plt.xlim(0,0.05e-7)
# plt.legend()

# # %% Prepare solver plane H
# solver = FDTD(input_filename = 'wireCoupling_long_H.fdtd.json', path_to_exe=SEMBA_EXE)
# # wood
# solver["materials"][5]["relativePermittivity"] = 4.0
# #styrofoam
# solver["materials"][6]["relativePermittivity"] = 1.1
# solver["materialAssociations"][1]["endTerminalId"] = 3
# solver.cleanUp()
# solver.run() 
# %% output, plane H
# plt.figure()
# plt.plot(v_inputH['time'], v_inputH['voltage_0'],'.-', label='input, H plane')
# plt.plot(v_outputH['time'], v_outputH['voltage_0'],'.-', label='output, H plane')
# plt.grid(which='both')
# plt.xlabel('time [s]')
# plt.ylabel('V [V]')
# plt.legend()

# # %% Prepare solver plane V
# solver = FDTD(input_filename = 'wireCoupling_long_V.fdtd.json', path_to_exe=SEMBA_EXE)
# # wood
# solver["materials"][5]["relativePermittivity"] = 4.0
# #styrofoam
# solver["materials"][6]["relativePermittivity"] = 1.1
# solver["materialAssociations"][1]["endTerminalId"] = 3
# solver.cleanUp()
# solver.run() 
# # %% output, plane v
# v_input_V =  Probe(solver.getSolvedProbeFilenames("circuit_probe_wire")[0])
# v_output_V =  Probe(solver.getSolvedProbeFilenames("circuit_probe_segment")[0])
# plt.figure()
# plt.plot(v_input_V['time'], v_input_V['voltage_0'],'.-', label='input, V plane')
# plt.plot(v_output_V['time'], v_output_V['voltage_0'],'.-', label='output, V plane')
# plt.grid(which='both')
# plt.xlabel('time [s]')
# plt.ylabel('V [V]')
# plt.legend()


# %% Postprocess

t = v_input['time']
dt = t[1] - t[0]
fq  = fftfreq(len(t))/dt
Vf_input = fft(v_input['voltage_0'])
Vf_input_H = fft(v_inputH['voltage_0'])
Vf_input_V = fft(v_inputV['voltage_0'])

Vf_output = fft(v_output['voltage_0'])
Vf_output_H = fft(v_outputH['voltage_0'])
Vf_output_V = fft(v_outputV['voltage_0'])

fmin = 80e6
fmax = 500e6
idx_min = (np.abs(fq - fmin)).argmin()
idx_max = (np.abs(fq - fmax)).argmin()
f = fq[idx_min:idx_max]
Vff_input = Vf_input[idx_min:idx_max]
Vff_output = Vf_output[idx_min:idx_max]
Vff_input_H = Vf_input_H[idx_min:idx_max]
Vff_output_H = Vf_output_H[idx_min:idx_max]
Vff_input_V = Vf_input_V[idx_min:idx_max]
Vff_output_V = Vf_output_V[idx_min:idx_max]

fig, ax = plt.subplots(1,2,figsize=(12, 4))
ax[0].plot(f*1e-6, 20*np.log10(np.abs(Vff_input_H)),'b.-', label='H plane, input')
ax[0].plot(f*1e-6, 20*np.log10(np.abs(Vff_input_V)),'r.-', label='V plane, input')
ax[0].plot(f*1e-6, 20*np.log10(np.abs(Vff_input)),'y.-', label='no plane, input')
ax[0].grid(which='both')
ax[0].set_xlabel('f [MHz]')
ax[0].set_ylabel('V [dB]')
ax[0].set_xlim(0.8e2,5e2)
# ax[0].set_ylim(40,100)
ax[0].legend()

ax[1].plot(f*1e-6, 20*np.log10(np.abs(Vff_output_H)),'b.-', label='H plane, output')
ax[1].plot(f*1e-6, 20*np.log10(np.abs(Vff_output_V)),'r.-', label='V plane, output')
ax[1].plot(f*1e-6, 20*np.log10(np.abs(Vff_output)),'y.-', label='no plane, output')
ax[1].grid(which='both')
ax[1].set_xlabel('f [MHz]')
ax[1].set_ylabel('V [dB]')
ax[1].set_xlim(0.8e2,5e2)
# ax[1].set_ylim(40,100)
ax[1].legend()

plt.show()

# %% Postprocess

n= 0 
for i in range(len(f)):
    if ((f[i]>0.8e8) and (f[i] < 5e8)):
        n = n+1
data_in_no_plane = np.zeros((n, 2))
data_in_H_plane = np.zeros((n, 2))
data_in_V_plane = np.zeros((n, 2))
data_out_no_plane = np.zeros((n, 2))
data_out_H_plane = np.zeros((n, 2))
data_out_V_plane = np.zeros((n, 2))
n = 0
for i in range(len(f)):
    if ((f[i]>0.8e8) and (f[i] < 5e8)):
        data_in_no_plane[n,0]  = f[i]
        data_in_H_plane[n,0]   = f[i]
        data_in_V_plane[n,0]   = f[i]
        data_out_no_plane[n,0] = f[i]
        data_out_H_plane[n,0]  = f[i]
        data_out_V_plane[n,0]  = f[i]

        data_in_no_plane[n,1] = 20*np.log10(np.abs(Vff_input[i]))
        data_in_H_plane[n,1]  = 20*np.log10(np.abs(Vff_input_H[i]))
        data_in_V_plane[n,1]  = 20*np.log10(np.abs(Vff_input_V[i]))
        data_out_no_plane[n,1] = 20*np.log10(np.abs(Vff_output[i]))
        data_out_H_plane[n,1] = 20*np.log10(np.abs(Vff_output_H[i]))
        data_out_V_plane[n,1] = 20*np.log10(np.abs(Vff_output_V[i]))
        n = n + 1
np.savetxt('input_no_plane_A.2.dat', data_in_no_plane)
np.savetxt('input_H_plane_A.2.dat', data_in_H_plane)
np.savetxt('input_V_plane_A.2.dat', data_in_V_plane)
np.savetxt('output_no_plane_A.2.dat', data_out_no_plane)
np.savetxt('output_H_plane_A.2.dat', data_out_H_plane)
np.savetxt('output_V_plane_A.2.dat', data_out_V_plane)


# %%
