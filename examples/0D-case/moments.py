#!/usr/bin/env python

import numpy as np
from os.path import join

dirname = "reference"
filename = join(dirname, "CSD-t0.dat")
L0, f0 = np.loadtxt(filename, unpack=True)


def calculate_moment(L, f, order=0):
    m = 0
    dL = np.diff(L)
    LMid = (L[1:] + L[:-1]) / 2
    fMid = (f[1:] + f[:-1]) / 2
    m = np.sum(fMid * LMid**order * dL)
    return m


rho = 1200
C = 0.6 * rho
kv = np.pi / 6
m3old = calculate_moment(L0, f0, 3)
G = 1e-8
t = 0
dt = 1
t_max = 2000
orders = range(6)

with open(join(dirname, "moments.dat"), "w") as f:
    f.write("t\tC")
    for order in orders:
        f.write(f"\tM{order}")
    f.write("\n")
    while t <= t_max:
        f.write(f"{t}\t{C}")
        L = L0 + G * t
        for order in orders:
            m = calculate_moment(L, f0, order)
            f.write(f"\t{m:e}")
            if order == 3:
                m3 = m
        f.write("\n")
        t += dt
        C -= kv * rho * (m3 - m3old)
        m3old = m3
