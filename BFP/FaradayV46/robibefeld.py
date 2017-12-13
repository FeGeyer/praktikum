import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

mhub = const.value('Bohr magneton')     # das gelibete Borhsche Magneton zeigt wie man Scipy Constants benutzt


def mittel(x):  # the real mean()-ing of life
    return ufloat(np.mean(x), np.std(x, ddof=1)/np.sqrt(len(x)))


def relf(l, m):  # in Prozent
    return (np.absolute(l-m)/l)*100


def fitf(x, a, b):
    return a*x**2 + b


z, B = np.genfromtxt('BFeld.txt', unpack=True)
# z = z-80    # Zentrum des Magfeldes
# z *= 10**-3
# B *= 10**-3

# Fit
params, cov = curve_fit(fitf, z, B)
errors = np.sqrt(np.diag(cov))
a = ufloat(params[0], errors[0])
b = ufloat(params[1], errors[1])
print("Hier a und b: ", a, b)
# Tabelle
# np.savetxt('BFeldtab.txt',np.column_stack([B,z]), delimiter=' & ',
#            newline= r'\\'+'\n' )

x = np.linspace(60, 100, 1000)
# plt.subplot(1, 2, 1)
plt.figure(1)
plt.plot(z, B, 'ro', label='Mag. Feldstärke')
plt.plot(x, fitf(x, *params), 'b-', label='Parabel')
plt.xlabel(r'$z \:/\: m$')
plt.ylabel(r'$B \:/\: T$')
plt.grid()
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('BFeld.pdf')
# plt.clf()
