#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reference: https://doi.org/10.1021/acs.jpcb.2c03897
"v8.3" was re-written in python script to run in university computing resources
There were some changes in loops.
jongcheol1422@gmail.com
"""

import numpy as np
import random
import pandas as pd
import os

#---------------------------------------------input1-----------------------------------------------
c = 299792458 #m/s
lambda_IR = 3500 * 10**-9 #m
lambda_Vis = 800 * 10**-9
lambda_SFG = (lambda_IR * lambda_Vis) / (lambda_IR + lambda_Vis)

n_water_IR = 1
n_water_Vis = 1
n_water_SFG = 1
n_Cellulose_IR = 1.4577
n_Cellulose_Vis = 1.464
n_Cellulose_SFG = 1.467

Iangle_Vis = 18
Iangle_IR = -18

#calc1-------------------------------------------------------------------------

Rangle_Vis = -Iangle_Vis
Rangle_IR = -Iangle_IR

Tangle_IR = np.arcsin((n_water_IR * np.sin(np.radians(Iangle_IR))) / n_Cellulose_IR) * (180 / np.pi) #Snell's law
Tangle_Vis = np.arcsin((n_water_Vis * np.sin(np.radians(Iangle_Vis))) / n_Cellulose_Vis) * (180 / np.pi)

#phase-matching condition along x-direction (momentum conservation)
Rangle_SFG = -np.arcsin(
    (n_water_Vis * (c / lambda_Vis) * np.sin(np.radians(Iangle_Vis)) +
     n_water_IR * (c / lambda_IR) * np.sin(np.radians(Iangle_IR))) / (n_water_SFG * (c / lambda_SFG))
) * (180 / np.pi)

Tangle_SFG = np.arcsin(
    (n_Cellulose_Vis * (c / lambda_Vis) * np.sin(np.radians(Tangle_Vis)) +
     n_Cellulose_IR * (c / lambda_IR) * np.sin(np.radians(Tangle_IR))) / (n_Cellulose_SFG * (c / lambda_SFG))
) * (180 / np.pi)

Iangle_SFG = -Rangle_SFG


#to check the medium (air/water) and propagation mode (co/counter)
if n_water_IR == 1:
    medium = "Air"
elif n_water_IR > 1:
    medium = "Water"
else:
    medium = "need to check"
    
if Iangle_IR > 0:
    mode = "Co-propagating"
elif Iangle_IR < 0:
    mode = "Counter-propagating"
else:
    mode = "need to check"

aaI = [Iangle_SFG, Iangle_Vis, Iangle_IR]
aaR = [Rangle_SFG, Rangle_Vis, Rangle_IR]
aaT = [Tangle_SFG, Tangle_Vis, Tangle_IR]    
aaI = np.array(aaI)
aaR = np.array(aaR)
aaT = np.array(aaT)

#Lc calculation - phase-matching condition along z-direction
kr1 = (2 * np.pi * n_water_IR / lambda_IR) * np.cos(np.radians(Iangle_IR))
kr2 = (2 * np.pi * n_water_Vis / lambda_Vis) * np.cos(np.radians(Iangle_Vis))
kr3 = (2 * np.pi * n_water_SFG / lambda_SFG) * np.cos(np.radians(Iangle_SFG))
dkr = kr3 + (kr1 + kr2)
Lcr = np.pi / dkr

kt1 = (2 * np.pi * n_Cellulose_IR / lambda_IR) * np.cos(np.radians(Tangle_IR))
kt2 = (2 * np.pi * n_Cellulose_Vis / lambda_Vis) * np.cos(np.radians(Tangle_Vis))
kt3 = (2 * np.pi * n_Cellulose_SFG / lambda_SFG) * np.cos(np.radians(Tangle_SFG))
dkt = kt3 - (kt1 + kt2)
Lct = np.pi / dkt

print("---------------------check laser/exp configuration----------------------")
print("Medium =", medium)
print("Propagation =", mode)
print("Incident Angle of SFG, Vis, and IR =", aaI)
print("Reflection Angle of SFG, Vis, and IR =", aaR)
print("Transmission Angle of SFG, Vis, and IR =", aaT)
print(f'Coherence Length (nm, Reflection) = {Lcr * 1e9:.0f} nm')
print(f'Coherence Length (um, Transmission) = {Lct * 1e6:.0f} um')


#calculation2-Lfactor (Fresnel coeff) or C_fres in Eqn-------------------------
RLxx = [
    (2 * n_water_SFG * np.cos(np.radians(Tangle_SFG))) /
    (n_water_SFG * np.cos(np.radians(Tangle_SFG)) + n_Cellulose_SFG * np.cos(np.radians(Iangle_SFG))),
    (2 * n_water_Vis * np.cos(np.radians(Tangle_Vis))) /
    (n_water_Vis * np.cos(np.radians(Tangle_Vis)) + n_Cellulose_Vis * np.cos(np.radians(Iangle_Vis))),
    (2 * n_water_IR * np.cos(np.radians(Tangle_IR))) /
    (n_water_IR * np.cos(np.radians(Tangle_IR)) + n_Cellulose_IR * np.cos(np.radians(Iangle_IR)))
]

RLyy = [
    (2 * n_water_SFG * np.cos(np.radians(Iangle_SFG))) /
    (n_water_SFG * np.cos(np.radians(Iangle_SFG)) + n_Cellulose_SFG * np.cos(np.radians(Tangle_SFG))),
    (2 * n_water_Vis * np.cos(np.radians(Iangle_Vis))) /
    (n_water_Vis * np.cos(np.radians(Iangle_Vis)) + n_Cellulose_Vis * np.cos(np.radians(Tangle_Vis))),
    (2 * n_water_IR * np.cos(np.radians(Iangle_IR))) /
    (n_water_IR * np.cos(np.radians(Iangle_IR)) + n_Cellulose_IR * np.cos(np.radians(Tangle_IR)))
]

RLzz = [
    ((2 * n_Cellulose_SFG * np.cos(np.radians(Iangle_SFG))) /
     (n_water_SFG * np.cos(np.radians(Tangle_SFG)) + n_Cellulose_SFG * np.cos(np.radians(Iangle_SFG))) *
     (n_water_SFG / n_Cellulose_SFG)**2),
    ((2 * n_Cellulose_Vis * np.cos(np.radians(Iangle_Vis))) /
     (n_water_Vis * np.cos(np.radians(Tangle_Vis)) + n_Cellulose_Vis * np.cos(np.radians(Iangle_Vis))) *
     (n_water_Vis / n_Cellulose_Vis)**2),
    ((2 * n_Cellulose_IR * np.cos(np.radians(Iangle_IR))) /
     (n_water_IR * np.cos(np.radians(Tangle_IR)) + n_Cellulose_IR * np.cos(np.radians(Iangle_IR))) *
     (n_water_IR / n_Cellulose_IR)**2)
]

TLxx = [
    (2 * n_Cellulose_SFG * np.cos(np.radians(Iangle_SFG))) /
    (n_water_SFG * np.cos(np.radians(Tangle_SFG)) + n_Cellulose_SFG * np.cos(np.radians(Iangle_SFG))),
    (2 * n_water_Vis * np.cos(np.radians(Tangle_Vis))) /
    (n_water_Vis * np.cos(np.radians(Tangle_Vis)) + n_Cellulose_Vis * np.cos(np.radians(Iangle_Vis))),
    (2 * n_water_IR * np.cos(np.radians(Tangle_IR))) /
    (n_water_IR * np.cos(np.radians(Tangle_IR)) + n_Cellulose_IR * np.cos(np.radians(Iangle_IR)))
]

TLyy = [
    (2 * n_Cellulose_SFG * np.cos(np.radians(Tangle_SFG))) /
    (n_water_SFG * np.cos(np.radians(Iangle_SFG)) + n_Cellulose_SFG * np.cos(np.radians(Tangle_SFG))),
    (2 * n_water_Vis * np.cos(np.radians(Iangle_Vis))) /
    (n_water_Vis * np.cos(np.radians(Iangle_Vis)) + n_Cellulose_Vis * np.cos(np.radians(Tangle_Vis))),
    (2 * n_water_IR * np.cos(np.radians(Iangle_IR))) /
    (n_water_IR * np.cos(np.radians(Iangle_IR)) + n_Cellulose_IR * np.cos(np.radians(Tangle_IR)))
]

TLzz = [
    ((2 * n_water_SFG * np.cos(np.radians(Tangle_SFG))) /
     (n_water_SFG * np.cos(np.radians(Tangle_SFG)) + n_Cellulose_SFG * np.cos(np.radians(Iangle_SFG))) *
     (n_Cellulose_SFG / n_Cellulose_SFG)**2),
    ((2 * n_Cellulose_Vis * np.cos(np.radians(Iangle_Vis))) /
     (n_water_Vis * np.cos(np.radians(Tangle_Vis)) + n_Cellulose_Vis * np.cos(np.radians(Iangle_Vis))) *
     (n_water_Vis / n_Cellulose_Vis)**2),
    ((2 * n_Cellulose_IR * np.cos(np.radians(Iangle_IR))) /
     (n_water_IR * np.cos(np.radians(Tangle_IR)) + n_Cellulose_IR * np.cos(np.radians(Iangle_IR))) *
     (n_water_IR / n_Cellulose_IR)**2)
]

RL = [RLxx, RLyy, RLzz]
TL = [TLxx, TLyy, TLzz]
RL = np.array(RL)
TL = np.array(TL)

#constant for projection
RIVSFG = [np.sign(Rangle_SFG) * np.cos(np.radians(Rangle_SFG)), 1, abs(np.sin(np.radians(Rangle_SFG)))]
RIVVIS = [np.cos(np.radians(Iangle_Vis)), 1, np.sin(np.radians(Iangle_Vis))]
RIVIR = [np.sign(Iangle_IR) * np.cos(np.radians(Iangle_IR)), 1, np.sign(Iangle_IR) * np.sin(np.radians(Iangle_IR))]

TIVSFG = [np.sign(Tangle_SFG) * np.cos(np.radians(Tangle_SFG)), 1, abs(np.sin(np.radians(Tangle_SFG)))]
TIVVIS = [np.cos(np.radians(Iangle_Vis)), 1, np.sin(np.radians(Iangle_Vis))]
TIVIR = [np.sign(Iangle_IR) * np.cos(np.radians(Iangle_IR)), 1, np.sign(Iangle_IR) * np.sin(np.radians(Iangle_IR))]

RIV = [RIVSFG, RIVVIS, RIVIR]
TIV = [TIVSFG, TIVVIS, TIVIR]
RIV = np.array(RIV)
TIV = np.array(TIV)

#calculate IR transition dipole vector
def calc_IRdipole(alpha, beta):
    return [np.sin(alpha) * np.cos(beta), np.sin(alpha) * np.sin(beta), np.cos(alpha)]

#dipole configurations
OHdipole = {'alpha': np.radians(30), 'beta': np.radians(0)}
CHdipole = {'alpha': np.radians(60), 'beta': np.radians(90)}

dipole_OH = calc_IRdipole(OHdipole['alpha'], OHdipole['beta'])
dipole_CH = calc_IRdipole(CHdipole['alpha'], CHdipole['beta'])

#hyperpolarizability tensors calculation
hyperOH = np.array([[[dipole_OH[k] * 0.9 for k in range(3)] for j in range(3)] for i in range(3)])
hyperCH = np.array([[[dipole_CH[k] * 1 for k in range(3)] for j in range(3)] for i in range(3)])

print("------------------1st check (copy&paste test1.txt)----------------------")



print("--------------------------1st check ends here---------------------------")

#eul
def euler_matrix(phi, theta, psi):  #euler matrix (ZYZ convention)
    phi = np.radians(phi)
    theta = np.radians(theta)
    psi = np.radians(psi)
    
    c_phi = np.cos(phi)
    s_phi = np.sin(phi)
    c_theta = np.cos(theta)
    s_theta = np.sin(theta)
    c_psi = np.cos(psi)
    s_psi = np.sin(psi)

    euler_mat = np.array([
        [c_theta * c_phi * c_psi - s_phi * s_psi, -c_psi * s_phi - c_theta * c_phi * s_psi, c_phi * s_theta],
        [c_theta * c_psi * s_phi + c_phi * s_psi, c_phi * c_psi - c_theta * s_phi * s_psi, s_phi * s_theta],
        [-c_psi * s_theta, s_theta * s_psi, c_theta]
    ])

    return euler_mat

#calculate Chai2 values (27 for each R and T).
#each of 27 hyper go through eul, Lfactor, and projection (IV) to result in one Chai2 value.
#eul, RL, RIV are all 3x3 matrices. 

def Rchai2OH(var1, var2, var3, phi, theta, psi):
    eul = euler_matrix(phi, theta, psi)
    return sum(
        eul[var1 - 1, i] * eul[var2 - 1, j] * eul[var3 - 1, k] * \
            RL[var1 - 1, 0] * RL[var2 - 1, 1] * RL[var3 - 1, 2] * \
        RIV[0, var1 - 1] * RIV[1, var2 - 1] * RIV[2, var3 - 1] * hyperOH[i, j, k]
        for i in range(3) for j in range(3) for k in range(3)
    )

def Rchai2CH(var1, var2, var3, phi, theta, psi):
    eul = euler_matrix(phi, theta, psi)
    return sum(
        eul[var1 - 1, i] * eul[var2 - 1, j] * eul[var3 - 1, k] * \
            RL[var1 - 1, 0] * RL[var2 - 1, 1] * RL[var3 - 1, 2] * \
        RIV[0, var1 - 1] * RIV[1, var2 - 1] * RIV[2, var3 - 1] * hyperCH[i, j, k]
        for i in range(3) for j in range(3) for k in range(3)
    )

def Tchai2OH(var1, var2, var3, phi, theta, psi):
    eul = euler_matrix(phi, theta, psi)
    return sum(
        eul[var1 - 1, i] * eul[var2 - 1, j] * eul[var3 - 1, k] * \
            TL[var1 - 1, 0] * TL[var2 - 1, 1] * TL[var3 - 1, 2] * \
        TIV[0, var1 - 1] * TIV[1, var2 - 1] * TIV[2, var3 - 1] * hyperOH[i, j, k]
        for i in range(3) for j in range(3) for k in range(3)
    )

def Tchai2CH(var1, var2, var3, phi, theta, psi):
    eul = euler_matrix(phi, theta, psi)
    return sum(
        eul[var1 - 1, i] * eul[var2 - 1, j] * eul[var3 - 1, k] * \
            TL[var1 - 1, 0] * TL[var2 - 1, 1] * TL[var3 - 1, 2] * \
        TIV[0, var1 - 1] * TIV[1, var2 - 1] * TIV[2, var3 - 1] * hyperCH[i, j, k]
        for i in range(3) for j in range(3) for k in range(3)
    )

#combination (sum) of chi2 values for each PC (polarizaiton combination) and dipole (CH&OH)
def compute_RSSPOH(phi, theta, psi):
    return (Rchai2OH(2, 2, 1, phi, theta, psi) + Rchai2OH(2, 2, 3, phi, theta, psi))

def compute_RSSPCH(phi, theta, psi):
    return (Rchai2CH(2, 2, 1, phi, theta, psi) + Rchai2CH(2, 2, 3, phi, theta, psi))

def compute_RPPPOH(phi, theta, psi):
    return (Rchai2OH(1, 1, 1, phi, theta, psi) + Rchai2OH(3, 1, 1, phi, theta, psi) + 
            Rchai2OH(1, 3, 1, phi, theta, psi) + Rchai2OH(1, 1, 3, phi, theta, psi) + 
            Rchai2OH(1, 3, 3, phi, theta, psi) + Rchai2OH(3, 1, 3, phi, theta, psi) + 
            Rchai2OH(3, 3, 1, phi, theta, psi) + Rchai2OH(3, 3, 3, phi, theta, psi))

def compute_RPPPCH(phi, theta, psi):
    return (Rchai2CH(1, 1, 1, phi, theta, psi) + Rchai2CH(3, 1, 1, phi, theta, psi) + 
            Rchai2CH(1, 3, 1, phi, theta, psi) + Rchai2CH(1, 1, 3, phi, theta, psi) + 
            Rchai2CH(1, 3, 3, phi, theta, psi) + Rchai2CH(3, 1, 3, phi, theta, psi) + 
            Rchai2CH(3, 3, 1, phi, theta, psi) + Rchai2CH(3, 3, 3, phi, theta, psi))

def compute_RPSPOH(phi, theta, psi):
    return (Rchai2OH(1, 2, 1, phi, theta, psi) + Rchai2OH(3, 2, 1, phi, theta, psi) + 
            Rchai2OH(1, 2, 3, phi, theta, psi) + Rchai2OH(3, 2, 3, phi, theta, psi))

def compute_RPSPCH(phi, theta, psi):
    return (Rchai2CH(1, 2, 1, phi, theta, psi) + Rchai2CH(3, 2, 1, phi, theta, psi) + 
            Rchai2CH(1, 2, 3, phi, theta, psi) + Rchai2CH(3, 2, 3, phi, theta, psi))

def compute_RSPPOH(phi, theta, psi):
    return (Rchai2OH(2, 1, 1, phi, theta, psi) + Rchai2OH(2, 3, 1, phi, theta, psi) + 
            Rchai2OH(2, 1, 3, phi, theta, psi) + Rchai2OH(2, 3, 3, phi, theta, psi))

def compute_RSPPCH(phi, theta, psi):
    return (Rchai2CH(2, 1, 1, phi, theta, psi) + Rchai2CH(2, 3, 1, phi, theta, psi) + 
            Rchai2CH(2, 1, 3, phi, theta, psi) + Rchai2CH(2, 3, 3, phi, theta, psi))

def compute_RSSSOH(phi, theta, psi):
    return Rchai2OH(2, 2, 2, phi, theta, psi)

def compute_RSSSCH(phi, theta, psi):
    return Rchai2CH(2, 2, 2, phi, theta, psi)

def compute_RPSSOH(phi, theta, psi):
    return (Rchai2OH(1, 2, 2, phi, theta, psi) + Rchai2OH(3, 2, 2, phi, theta, psi))

def compute_RPSSCH(phi, theta, psi):
    return (Rchai2CH(1, 2, 2, phi, theta, psi) + Rchai2CH(3, 2, 2, phi, theta, psi))

def compute_RPPSOH(phi, theta, psi):
    return (Rchai2OH(1, 1, 2, phi, theta, psi) + Rchai2OH(3, 1, 2, phi, theta, psi) + 
            Rchai2OH(1, 3, 2, phi, theta, psi) + Rchai2OH(3, 3, 2, phi, theta, psi))

def compute_RPPSCH(phi, theta, psi):
    return (Rchai2CH(1, 1, 2, phi, theta, psi) + Rchai2CH(3, 1, 2, phi, theta, psi) + 
            Rchai2CH(1, 3, 2, phi, theta, psi) + Rchai2CH(3, 3, 2, phi, theta, psi))

def compute_RSPSOH(phi, theta, psi):
    return (Rchai2OH(2, 1, 2, phi, theta, psi) + Rchai2OH(2, 3, 2, phi, theta, psi))

def compute_RSPSCH(phi, theta, psi):
    return (Rchai2CH(2, 1, 2, phi, theta, psi) + Rchai2CH(2, 3, 2, phi, theta, psi))

def compute_TSSPOH(phi, theta, psi):
    return (Tchai2OH(2, 2, 1, phi, theta, psi) + Tchai2OH(2, 2, 3, phi, theta, psi))

def compute_TSSPCH(phi, theta, psi):
    return (Tchai2CH(2, 2, 1, phi, theta, psi) + Tchai2CH(2, 2, 3, phi, theta, psi))

def compute_TPPPOH(phi, theta, psi):
    return (Tchai2OH(1, 1, 1, phi, theta, psi) + Tchai2OH(3, 1, 1, phi, theta, psi) + 
            Tchai2OH(1, 3, 1, phi, theta, psi) + Tchai2OH(1, 1, 3, phi, theta, psi) + 
            Tchai2OH(1, 3, 3, phi, theta, psi) + Tchai2OH(3, 1, 3, phi, theta, psi) + 
            Tchai2OH(3, 3, 1, phi, theta, psi) + Tchai2OH(3, 3, 3, phi, theta, psi))

def compute_TPPPCH(phi, theta, psi):
    return (Tchai2CH(1, 1, 1, phi, theta, psi) + Tchai2CH(3, 1, 1, phi, theta, psi) + 
            Tchai2CH(1, 3, 1, phi, theta, psi) + Tchai2CH(1, 1, 3, phi, theta, psi) + 
            Tchai2CH(1, 3, 3, phi, theta, psi) + Tchai2CH(3, 1, 3, phi, theta, psi) + 
            Tchai2CH(3, 3, 1, phi, theta, psi) + Tchai2CH(3, 3, 3, phi, theta, psi))

def compute_TPSPOH(phi, theta, psi):
    return (Tchai2OH(1, 2, 1, phi, theta, psi) + Tchai2OH(3, 2, 1, phi, theta, psi) + 
            Tchai2OH(1, 2, 3, phi, theta, psi) + Tchai2OH(3, 2, 3, phi, theta, psi))

def compute_TPSPCH(phi, theta, psi):
    return (Tchai2CH(1, 2, 1, phi, theta, psi) + Tchai2CH(3, 2, 1, phi, theta, psi) + 
            Tchai2CH(1, 2, 3, phi, theta, psi) + Tchai2CH(3, 2, 3, phi, theta, psi))

def compute_TSPPOH(phi, theta, psi):
    return (Tchai2OH(2, 1, 1, phi, theta, psi) + Tchai2OH(2, 3, 1, phi, theta, psi) + 
            Tchai2OH(2, 1, 3, phi, theta, psi) + Tchai2OH(2, 3, 3, phi, theta, psi))

def compute_TSPPCH(phi, theta, psi):
    return (Tchai2CH(2, 1, 1, phi, theta, psi) + Tchai2CH(2, 3, 1, phi, theta, psi) + 
            Tchai2CH(2, 1, 3, phi, theta, psi) + Tchai2CH(2, 3, 3, phi, theta, psi))

def compute_TSSSOH(phi, theta, psi):
    return Tchai2OH(2, 2, 2, phi, theta, psi)

def compute_TSSSCH(phi, theta, psi):
    return Tchai2CH(2, 2, 2, phi, theta, psi)

def compute_TPSSOH(phi, theta, psi):
    return (Tchai2OH(1, 2, 2, phi, theta, psi) + Tchai2OH(3, 2, 2, phi, theta, psi))

def compute_TPSSCH(phi, theta, psi):
    return (Tchai2CH(1, 2, 2, phi, theta, psi) + Tchai2CH(3, 2, 2, phi, theta, psi))

def compute_TPPSOH(phi, theta, psi):
    return (Tchai2OH(1, 1, 2, phi, theta, psi) + Tchai2OH(3, 1, 2, phi, theta, psi) + 
            Tchai2OH(1, 3, 2, phi, theta, psi) + Tchai2OH(3, 3, 2, phi, theta, psi))

def compute_TPPSCH(phi, theta, psi):
    return (Tchai2CH(1, 1, 2, phi, theta, psi) + Tchai2CH(3, 1, 2, phi, theta, psi) + 
            Tchai2CH(1, 3, 2, phi, theta, psi) + Tchai2CH(3, 3, 2, phi, theta, psi))

def compute_TSPSOH(phi, theta, psi):
    return (Tchai2OH(2, 1, 2, phi, theta, psi) + Tchai2OH(2, 3, 2, phi, theta, psi))

def compute_TSPSCH(phi, theta, psi):
    return (Tchai2CH(2, 1, 2, phi, theta, psi) + Tchai2CH(2, 3, 2, phi, theta, psi))

compute_functions = {
    "RSSPOH": lambda phi, theta, psi: compute_RSSPOH(phi, theta, psi),
    "RSSPCH": lambda phi, theta, psi: compute_RSSPCH(phi, theta, psi),
    "RPPPOH": lambda phi, theta, psi: compute_RPPPOH(phi, theta, psi),
    "RPPPCH": lambda phi, theta, psi: compute_RPPPCH(phi, theta, psi),
    "RPSPOH": lambda phi, theta, psi: compute_RPSPOH(phi, theta, psi),
    "RPSPCH": lambda phi, theta, psi: compute_RPSPCH(phi, theta, psi),
    "RSPPOH": lambda phi, theta, psi: compute_RSPPOH(phi, theta, psi),
    "RSPPCH": lambda phi, theta, psi: compute_RSPPCH(phi, theta, psi),
    "RSSSOH": lambda phi, theta, psi: compute_RSSSOH(phi, theta, psi),
    "RSSSCH": lambda phi, theta, psi: compute_RSSSCH(phi, theta, psi),
    "RPSSOH": lambda phi, theta, psi: compute_RPSSOH(phi, theta, psi),
    "RPSSCH": lambda phi, theta, psi: compute_RPSSCH(phi, theta, psi),
    "RPPSOH": lambda phi, theta, psi: compute_RPPSOH(phi, theta, psi),
    "RPPSCH": lambda phi, theta, psi: compute_RPPSCH(phi, theta, psi),
    "RSPSOH": lambda phi, theta, psi: compute_RSPSOH(phi, theta, psi),
    "RSPSCH": lambda phi, theta, psi: compute_RSPSCH(phi, theta, psi),
    "TSSPOH": lambda phi, theta, psi: compute_TSSPOH(phi, theta, psi),
    "TSSPCH": lambda phi, theta, psi: compute_TSSPCH(phi, theta, psi),
    "TPPPOH": lambda phi, theta, psi: compute_TPPPOH(phi, theta, psi),
    "TPPPCH": lambda phi, theta, psi: compute_TPPPCH(phi, theta, psi),
    "TPSPOH": lambda phi, theta, psi: compute_TPSPOH(phi, theta, psi),
    "TPSPCH": lambda phi, theta, psi: compute_TPSPCH(phi, theta, psi),
    "TSPPOH": lambda phi, theta, psi: compute_TSPPOH(phi, theta, psi),
    "TSPPCH": lambda phi, theta, psi: compute_TSPPCH(phi, theta, psi),
    "TSSSOH": lambda phi, theta, psi: compute_TSSSOH(phi, theta, psi),
    "TSSSCH": lambda phi, theta, psi: compute_TSSSCH(phi, theta, psi),
    "TPSSOH": lambda phi, theta, psi: compute_TPSSOH(phi, theta, psi),
    "TPSSCH": lambda phi, theta, psi: compute_TPSSCH(phi, theta, psi),
    "TPPSOH": lambda phi, theta, psi: compute_TPPSOH(phi, theta, psi),
    "TPPSCH": lambda phi, theta, psi: compute_TPPSCH(phi, theta, psi),
    "TSPSOH": lambda phi, theta, psi: compute_TSPSOH(phi, theta, psi),
    "TSPSCH": lambda phi, theta, psi: compute_TSPSCH(phi, theta, psi)
}

def gen_ang(phi1, theta1, psi, sDphi, sDtheta, sDpsi):
    if psi == 'random':
        gen_psi = float(random.uniform(0, 360))
    else:
        gen_psi = float(random.gauss(psi, sDpsi))
    
    return {
        'phi': float(random.gauss(phi1, sDphi)),
        'theta': float(random.gauss(theta1, sDtheta)),
        'psi': gen_psi
    }

print("------------------2nd check (copy&paste test2.txt)----------------------")



print("--------------------------2nd check ends here---------------------------")

#---------------------------------------------input2-----------------------------------------------
mode = 2 #1 for reflection 2 for transmission
dE = [0]
theta_angles = [90]
sDphi = [7.5]
sDtheta = [7.5]
dzlist = [5]
ulist = [4]

iteration = 1000
m = 12
#--------------------------------------------------------------------------------------------------

### for tests
tt1_a = np.zeros((25, iteration, m), dtype=np.complex128) # m OH dipoles, 5 iteration, at 25 phi angles
tt1_b = np.zeros((25, iteration), dtype=np.complex128)
tt1_c = np.zeros((25, iteration), dtype=np.complex128)
tt1_d = np.zeros(25,dtype=np.complex128)
#Ori angles
tt1_e = np.zeros((25, iteration, m, 3), dtype=np.float64) # m OH dipoles, 5 iteration, at 25 phi angles
#seq
tt1_f = np.zeros((25, iteration, m), dtype=np.float64) # m OH dipoles, 5 iteration, at 25 phi angles
###

pc = 2
while pc <= 2:
    if mode == 1:
        polarization_combination = ["RSSP", "RPPS", "RPPP", "RSSS", "RPSP", "RSPS", "RSPP", "RPSS"]
    else:
        polarization_combination = ["TSSP", "TPPS", "TPPP", "TSSS", "TPSP", "TSPS", "TSPP", "TPSS"]

    pcOH = polarization_combination[pc - 1] + "OH"
    pcCH = polarization_combination[pc - 1] + "CH"

    dk = dkr if mode == 1 else dkt

    dEN = 1
    dEList = dE
    while dEN <= len(dEList):
        theta_cases = theta_angles
        theta_casesN = 1
        while theta_casesN <= len(theta_cases):
            theta_fix = theta_cases[theta_casesN - 1]
            dz_cases = dzlist
            dz_casesN = 1
            while dz_casesN <= len(dz_cases):
                dz = dzlist[dz_casesN - 1]
                u_cases = ulist
                u_casesN = 1
                while u_casesN <= len(u_cases):
                    u = ulist[u_casesN - 1]
                    sDtheta_list = sDtheta
                    sDthetaN = 1
                    while sDthetaN <= len(sDtheta_list):
                        sDphi_list = sDphi
                        sDphiN = 1
                        while sDphiN <= len(sDphi_list):
                            phi_temp = 0
                            temp3 = []                            
                            while phi_temp <= 360:
                                phi1 = phi_temp
                                theta1 = theta_fix
                                ###
                                phi_index = int(phi_temp / 15) ### for tests
                                ###
                                
                                seqtemp = np.concatenate([
                                    np.ones(int(np.ceil(-0.1 + m * (dEList[dEN - 1] + 1) / 2))),
                                    np.zeros(int(np.floor(0.1 + m * (1 - (dEList[dEN - 1] + 1) / 2))))
                                ])
                                ###
                                #print(seqtemp)

                                G = np.zeros(iteration)
                                J = np.zeros(iteration)
                                ddd = np.zeros((111, 3))
                                fff = np.zeros((111, 3))

                                p = 1
                                q = 1

                                ddd[(p - 1) * 111 + q - 1, 1] = u
                                ddd[(p - 1) * 111 + q - 1, 0] = dz
                                fff[(p - 1) * 111 + q - 1, 1] = u
                                fff[(p - 1) * 111 + q - 1, 0] = dz

                                for w in range(1, iteration + 1):
                                    seq = random.sample(list(seqtemp), len(seqtemp))
                                    seq = np.array(seq)

                                    OHtemp1 = np.zeros(m, dtype=np.complex128)
                                    CHtemp1 = np.zeros(m, dtype=np.complex128)                                    
                                    for i in range(m):
                                        Ori1 = gen_ang(phi_temp, theta_fix, 'random', sDphi_list[sDphiN-1], sDtheta_list[sDthetaN-1], 90) #90 is random number
                                        Oritemp1 = Ori1
                                        Oritemp2 = {
                                            'phi': Oritemp1['phi'],
                                            'theta': Oritemp1['theta'] + 180,
                                            'psi': Oritemp1['psi'] + 180
                                        }
                                        ###
                                        #print(Oritemp1)
                                        
                                        Oritemp3 = Oritemp1 if seq[i] == 1 else Oritemp2

                                        phi_calc = float(Oritemp3['phi'])
                                        theta_calc = float(Oritemp3['theta'])
                                        psi_calc = float(Oritemp3['psi'])
                                        ###
                                        tt1_e[phi_index, w-1, i, 0] = phi_calc
                                        tt1_e[phi_index, w-1, i, 1] = theta_calc
                                        tt1_e[phi_index, w-1, i, 2] = psi_calc
                                        ###
                                        
                                        OHtemp2 = ((np.exp(-1j * dk * u * 1e-9) - 1) / dk) * \
                                            compute_functions[pcOH](phi_calc, theta_calc, psi_calc) * \
                                            np.exp(-1j * dk * (i * u * 1e-9 + dz * 1e-9))

                                        CHtemp2 = ((np.exp(-1j * dk * u * 1e-9) - 1) / dk) * \
                                            compute_functions[pcCH](phi_calc, theta_calc, psi_calc) * \
                                            np.exp(-1j * dk * (i * u * 1e-9 + dz * 1e-9))

                                        OHtemp1[i] = OHtemp2
                                        CHtemp1[i] = CHtemp2
                                        ###

                                        #--------For loop for i (m=12) ends here-----------------------------------------------

                                    tempOH = np.sum(OHtemp1)
                                    tempCH = np.sum(CHtemp1)
                                    G[w-1] = (np.conj(tempOH) * tempOH) / (np.cos(np.radians(Tangle_SFG)) ** 2) # G #elements=iteration
                                    J[w-1] = (np.conj(tempCH) * tempCH) / (np.cos(np.radians(Tangle_SFG)) ** 2)
                                    ###
                                    tt1_a[phi_index, w-1, :] = OHtemp1 # 12 OH dipoles at 5 iteration at 25 degrees
                                    tt1_b[phi_index, w-1] = tempOH # summed OH dipoles at 5 iteration at 25 degrees
                                    tt1_c[phi_index, w-1] = G[w-1] # 5 calculated OH intensity at each iteration at 25 degrees
                                    ###              
                                    tt1_f[phi_index, w-1] = seq # 12 OH dipoles at 5 iteration at 25 degrees
                                    #-------------For loop for w (iteration) ends here----------------------------------------
                                    
                                ddd[(p-1) * 111 + q-1, 2] = np.mean(G) # the # of elements in G: iteration
                                fff[(p-1) * 111 + q-1, 2] = np.mean(J)
                                ###
                                tt1_d[phi_index] = np.mean(G)
                                ###
                                
                                OHCHpara = [
                                    phi_temp,  
                                    theta_fix,
                                    psi_calc,  
                                    ddd[p-1, 0],
                                    ddd[p-1, 1],
                                    ddd[p-1, 2] / fff[p-1, 2] if fff[p-1, 2] > 0 else 0,
                                    ddd[p-1, 2],
                                    fff[p-1, 2]
                                ]

                                temp3.append(OHCHpara)  

                                phi_temp += 15                        
                                #---------------While loop for 0 ~360 phi ration ends here------------------------------------

                            filename = (f"DE{int(dEList[dEN - 1] * 100)}"
                                        f"-{polarization_combination[pc - 1]}"
                                        f"-Azi-sDphi{sDphi_list[sDphiN - 1]}"
                                        f"-sDth{sDtheta_list[sDthetaN - 1]}"
                                        f"-(azi{phi1:.0f}+tilt{theta_fix:.0f}+psi{'Q'})"
                                        f"-(m{round(m, 1)}+u{round(u, 1)}"
                                        f"+dz{round(dz, 1)}+int{iteration}).csv")

                            df = pd.DataFrame(temp3)
                            export_path = os.path.join(os.getcwd(), filename)
                            df.to_csv(export_path, index=False)

                            sDphiN += 1
                        sDthetaN += 1
                    u_casesN += 1
                dz_casesN += 1
            theta_casesN += 1
        dEN += 1
    pc += 1

print("------------------3rd check (copy&paste test3.txt)----------------------")

print("=======Angles=======")
""""
tt1_e (Ori1):  At each phi_temp angle (0-24), at each iteration (0-999), all 12 OH dipoles (0-11)"""
print(tt1_e)
print(tt1_e[24])
print(tt1_e.shape)

print("=======seq=======")
"""
"tt1_f (seq):  At each phi_temp angle (0-24), at each iteration (0-999), all 12 dipole directions (0-11)"""
print(tt1_f)
print(tt1_f[24])
print(tt1_f.shape)

print("=======Dipoles, int, mean =======")
"""
tt1_a (OHtemp1)  :  At each phi_temp angle (0-24), at each iteration (0-999), all 12 OH dipoles (0-11)
tt1_b (tempOH)   :  At each phi_temp angle (0-24), at each iteration (0-999), summed OH diples (.)
tt1_c (G[w-1])   :  At each phi_temp angle (0-24), at each iteration (0-999), calculated int (.)
tt1_d (np.mean.G):  At each phi_temp angle (0-24), mean intensity (.)
"""
print(tt1_c[24])
print(tt1_c.shape)
print(tt1_d)
print(tt1_d.shape)

print("--------------------------3rd check ends here---------------------------")

print("calculation completed")