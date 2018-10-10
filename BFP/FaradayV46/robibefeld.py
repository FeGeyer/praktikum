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

def mittel(x):  # the real mean()-ing of life
    return ufloat(np.mean(x), np.std(x, ddof=1)/np.sqrt(len(x)))


def relf(l, m):  # in Prozent
    return (np.absolute(l-m)/l)*100


def fitf(x, a, b, c):
    return a*x**2 + b*x + c


z, zrel, B = np.genfromtxt('BFeld.txt', unpack=True)
# z *= 10**-3
# B *= 10**-3

# Fit
params, cov = curve_fit(fitf, zrel, B)
errors = np.sqrt(np.diag(cov))
a = ufloat(params[0], errors[0])
b = ufloat(params[1], errors[1])
c = ufloat(params[2], errors[2])
print("Hier a, b und c: ", a, b, c)

# Tabelle
np.savetxt('TexTabellen/BFeldtab.txt',np.column_stack([B, z, zrel]), delimiter=' & ',
           newline= r' \\'+'\n', fmt='%.0f')

x = np.linspace(-30, 30, 1000)
y = fitf(x, a, b, c)
print("Maximum des B-Feldes: ", np.max(y))
# plt.subplot(1, 2, 1)
plt.figure(1)
plt.plot(zrel, B, 'rx', label=r'$B(z)$')
plt.plot(x, fitf(x, *params), 'b-', label='Parabel')
plt.xlabel(r'$z \:/\: \mathrm{mm}$')
plt.ylabel(r'$B \:/\: \mathrm{T}$')
plt.xlim(-20, 20)
plt.ylim(0, B.max()+25)
plt.grid()
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('BFeld.pdf')
# plt.clf()
