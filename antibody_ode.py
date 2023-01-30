# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 22:24:23 2023

@author: stefa
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# function that returns dz/dt
def model(concs,t,k):
    
    k1, k2, k3, k4, k5, k6, k7, k8 = k
    cI, rPoZ, P_cI, cI_dna, A, A_dna = concs
    
    dcI_dt = -k1*cI*rPoZ + k2*A -k3*cI*P_cI + k4*cI_dna
    
    drpoz_dt = -k1*cI*rPoZ + k2*A - k5*cI_dna*rPoZ + k6*A_dna
    
    dcIdna_dt = k3*cI*P_cI - k5*cI_dna*rPoZ + k6*A_dna
    
    dPcI_dt = -k3*cI*P_cI + k4*cI_dna - k7*A_dna*P_cI
    
    dA_dt = k1*cI*rPoZ - k2*A - k7*A_dna*P_cI + k8*A_dna
    
    dAdna_dt = k7*A*P_cI + k5*cI_dna - k8*A_dna
    
    dzdt = [dcI_dt, drpoz_dt, dPcI_dt, dcIdna_dt, dA_dt, dAdna_dt]
    return dzdt

# initial condition
concs = [10, 1, 0.1, 0, 0, 0]
# [cI, rPoZ, P_cI, cI_dna, A, A_dna]

k1 = 5 # rPoZ + cI binding
k2 = 1 # rPoZ + cI disocciation

k3 = 1 # cI + lambda repressor DNA binding
k4 = 1 # cI + lambda repressor DNA disocciation

k5 = k1 # rPoZ + cI binding w/ cI on DNA
k6 = k2 # rPoZ + cI disocciation w/ cI on DNA

k7 = k3 # cI + lambda repressor DNA binding w/ rPoZ on cI
k8 = k4 # cI + lambda repressor DNA disocciation w/ rPoZ on cI


k = [k1, k2, k3, k4, k5, k6, k7, k8]

# time points
t = np.linspace(0, 60, num = 10000)

# solve ODE
z, o = odeint(model,concs,t,args=(k,), full_output = 1)

# [cI, rPoZ, P_cI, cI_dna, A, A_dna]
plt.plot(t,z[:,5])

