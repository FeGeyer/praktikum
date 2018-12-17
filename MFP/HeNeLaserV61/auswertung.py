import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sem
# import uncertainties.unumpy as unp
# import pandas as pd
# from scipy.optimize import curve_fit
from uncertainties import ufloat
from scipy import constants

# Use latex fonts and text
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# pre-assignment


def bed1(L, r1, r2):
    return ((1 - L / r1) * (1 - L / r2))


def bed2(L, r1):
    return (1 - L / r1)


r1_1 = 1.4
r2_1 = 1.4

r1_2 = 1.4

x = np.linspace(0, 3)
a = bed1(x, r1_1, r2_1)
b = bed2(x, r1_2)

plt.plot(x, bed1(x, r1_1, r2_1), 'r-', label=r"$r_1$ = 0,5 m, $r_2$ = 0,2 m")
plt.plot(x, bed2(x, r1_2), 'b-', label=r"$r_1$ = 1 m, $r_2$ = 1,5 m")
plt.axhline(y=1, color='k', linestyle='--')
plt.axhline(y=0, color='k', linestyle='--')
plt.xlabel(r'Resonatorlaenge $L$ / m')
plt.ylabel(r'Stabilitaetsparameter $g_1 \cdot g_2$')
plt.ylim(-1, 1.7)
plt.xlim(0, 2.8)
plt.legend(loc="best")
plt.tight_layout
plt.savefig('Auswertung/Plots/g1g2.pdf')

print('----------------------------------------------------------------------')

# Wert in m
l1 = 0.6835
l2 = 1.295
l3 = 1.375

# Werte in MHz
P1 = np.array([229, 462, 691, 924, 1153, 1382])
P2 = np.array([118, 233, 351, 470, 584, 703, 817, 936, 1054, 1169, 1287])
P3 = np.array([110, 221, 331, 442, 553, 663, 774, 884, 995, 1102, 1212])

P1 = P1 * 10**6
P2 = P2 * 10**6
P3 = P3 * 10**6

nu1 = np.empty(len(P1) - 1)
nu2 = np.empty(len(P2) - 1)
nu3 = np.empty(len(P3) - 1)

for i in range(len(P1) - 1):
    nu1[i] = P1[i + 1] - P1[i]
for i in range(len(P2) - 1):
    nu2[i] = P2[i + 1] - P2[i]
for i in range(len(P3) - 1):
    nu3[i] = P3[i + 1] - P3[i]

dnu1 = ufloat(np.mean(nu1), sem(nu1))
dnu2 = ufloat(np.mean(nu2), sem(nu2))
dnu3 = ufloat(np.mean(nu3), sem(nu3))

# Some constants and values in SI units
c = constants.c
k = constants.k
T = 300
M = 20.1797 * constants.physical_constants["atomic mass constant"][0]
lam = 632.8 * 10**(-9)
f = c / lam

# v_mean from boltzman formular
v_mean = np.sqrt(2 * k * T / M)

# delta factor from doppler-effect
delta = (c + v_mean) / (c - v_mean)

# print frequency from 632.8nm wavelength and resulting doppler-effect
print(f, np.abs(2 * (f - delta * f)))

# print mean frequency between modes
print(dnu1, dnu2, dnu3)

print('----------------------------------------------------------------------')

print('Alles ausgeführt!')
