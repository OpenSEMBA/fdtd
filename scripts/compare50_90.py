#!/usr/bin/env python3
import subprocess

bin = '/home/luis/ugrfdtd/publico/cpp_build_nomtln/bin/cpp_tests'
for filt in [
    'PlanewavePwInBox.ShortRunProbeParity_First50Steps',
    'PlanewavePwInBox.MediumRunProbeParity_First90Steps',
]:
    subprocess.run([bin, f'--gtest_filter={filt}'], cwd='/home/luis/ugrfdtd/publico')

def read(p):
    f=[]
    for line in open(p):
        x=line.split()
        if len(x)>=2:
            try: f.append(float(x[1]))
            except: pass
    return f

a50 = read('/tmp/semba_pw_in_box_test/pw-in-box.fdtd_after_Ex_3_3_5.dat')
subprocess.run([bin, '--gtest_filter=PlanewavePwInBox.MediumRunProbeParity_First90Steps'], cwd='/home/luis/ugrfdtd/publico')
a90 = read('/tmp/semba_pw_in_box_test/pw-in-box.fdtd_after_Ex_3_3_5.dat')
print('len50',len(a50),'len90',len(a90))
for s in [50,51,90]:
    if s < len(a50) and s < len(a90):
        print(s, a50[s], a90[s], 'diff', abs(a50[s]-a90[s]))
