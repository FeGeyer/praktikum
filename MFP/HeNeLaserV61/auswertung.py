import numpy as np
import matplotlib.pyplot as plt
# import uncertainties.unumpy as unp
# import pandas as pd
# from scipy.optimize import curve_fit
# from uncertainties import ufloat
# from scipy import constants

# Use latex fonts and text
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# pre-assignment


def g1g2(L, r1, r2):
    return ((1 - L/r1)*(1 - L/r2))


r1_1 = 0.5
r2_1 = 0.2

r1_2 = 1
r2_2 = 1.5

x = np.linspace(0, 2.5)
a = g1g2(x, r1_1, r2_1)
a_cut = x[a >= 1][1]
b = g1g2(x, r1_2, r2_2)
b_cut = x[b >= 1][1]

print(a_cut)
print(b_cut)


plt.plot(x, g1g2(x, r1_1, r2_1), 'r.', label=r"$r_1$ = 0,5 m, $r_2$ = 0,2 m")
plt.plot(x, g1g2(x, r1_2, r2_2), 'b.', label=r"$r_1$ = 1 m, $r_2$ = 1,5 m")
plt.axhline(y=1, color='k', linestyle='--')
plt.axhline(y=0, color='k', linestyle='--')
plt.xlabel(r'Resonatorlaenge $L$ / m')
plt.ylabel(r'Stabilitaetsparameter $g_1 \cdot g_2$')
plt.ylim(-0.5, 1.7)
plt.legend(loc="best")
plt.tight_layout
plt.savefig('Auswertung/Plots/g1g2.pdf')

print('Alles ausgeführt!')
