#!/usr/bin/env python3
import subprocess
import sys

def read(path):
    f, inc = [], []
    for line in open(path):
        p = line.split()
        if len(p) < 3:
            continue
        try:
            f.append(float(p[1]))
            inc.append(float(p[2]))
        except ValueError:
            pass
    return f, inc

def corr(a, b):
    import numpy as np
    a = np.array(a) - np.mean(a)
    b = np.array(b) - np.mean(b)
    da = np.dot(a, a)
    db = np.dot(b, b)
    return 0.0 if da <= 0 or db <= 0 else float(np.dot(a, b) / np.sqrt(da * db))

bin_path = sys.argv[1] if len(sys.argv) > 1 else './cpp_build_nomtln/bin/cpp_tests'
subprocess.run([
    bin_path,
    '--gtest_filter=PlanewavePwInBox.MediumRunProbeParity_First90Steps:PlanewavePwInBox.Step100ProbeParity',
], cwd='/home/luis/ugrfdtd/publico', check=False)

ref = '/home/luis/ugrfdtd/publico/testData/cases/planewave'
got = '/tmp/semba_pw_in_box_test'
fb, _ = read(f'{got}/pw-in-box.fdtd_before_Ex_3_3_1.dat')
fi, ii = read(f'{got}/pw-in-box.fdtd_inbox_Ex_3_3_3.dat')
fa, _ = read(f'{got}/pw-in-box.fdtd_after_Ex_3_3_5.dat')
rfb, _ = read(f'{ref}/pw-in-box.fdtd_before_Ex_3_3_1.dat')
rfa, _ = read(f'{ref}/pw-in-box.fdtd_after_Ex_3_3_5.dat')
n = len(fa)
print('samples', n, 'corr', corr(fi, ii))
for s in [89, 90, 100]:
    print(f's{s}: b={fb[s]:.4e}/{rfb[s]:.4e} a={fa[s]:.4e}/{rfa[s]:.4e}')
err = 0
if corr(fi, ii) <= 0.999:
    err += 1
if any(abs(x) > 5e-4 for x in fb):
    err += 2
if any(abs(x) > 5e-4 for x in fa):
    err += 4
print('err', err, 'max_b', max(map(abs, fb)), 'max_a', max(map(abs, fa)))
