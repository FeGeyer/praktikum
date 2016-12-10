import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 10

#Einlesen der Daten
Phase, Uunv, Uver = np.genfromtxt('dataLockIn.txt', unpack=True)
AbsCm, ULED = np.genfromtxt('dataAbs.txt', unpack=True)
Gain2u3 = 5
Gain4 = 200

#Fitten von 3 und 4
def f(phi, A, dphi):
    return A * np.cos(phi+dphi)
params, covariance = curve_fit(f, (Phase / 180)*np.pi, Uunv/Gain2u3)
errors = np.sqrt(np.diag(covariance))
Anull = ufloat(params[0], errors[0])
deltaphi = ufloat(params[1], errors[1])

print(Anull)
print(deltaphi)
np.savetxt('Teil1.txt', np.column_stack([(Phase / 180), Uunv, Uver]),fmt="%.2f")
#plot
x_plot = np.linspace(0, 2*np.pi)

plt.figure(0)
plt.title("Auftragen von U gegen $\phi$ und Regression")

plt.ylabel('$U/V$')
plt.xlabel("$\phi/rad$")
plt.plot((Phase / 180)*np.pi, Uunv/Gain2u3, 'r+', label="Messwerte")
plt.plot(x_plot, f(x_plot, *params), 'b-', label='Regression')
plt.xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi],
           [r"$0$", r"$\frac{1}{2}\pi$", r"$\pi$", r"$\frac{3}{2}\pi$", r"$2\pi$"])
plt.legend(loc="best")
plt.tight_layout
plt.savefig('UvonPhi.pdf')

plt.figure(1)
plt.title("Spannung als Funktion des Abstandes zwischen LED und PD, doppelt logarithmisch")
plt.yscale("log")
plt.xscale("log")
plt.plot(AbsCm, ULED, 'b+', label='$U(r)$')
plt.legend(loc="best")
plt.tight_layout
plt.savefig('Abst.pdf')
