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
AbsCmWeg, ULEDweg = np.genfromtxt('dataWeg.txt', unpack=True)
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
#plt.title("Auftragen von U gegen $\phi$ und Regression")

plt.ylabel('$U/V$')
plt.xlabel("$\phi/rad$")
plt.plot((Phase / 180)*np.pi, Uunv/Gain2u3, 'r+', label="Messwerte")
plt.plot(x_plot, f(x_plot, *params), 'b-', label='Regression')
plt.xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi],
           [r"$0$", r"$\frac{1}{2}\pi$", r"$\pi$", r"$\frac{3}{2}\pi$", r"$2\pi$"])
plt.xlim(0, 2*np.pi)
plt.legend(loc="best")
plt.tight_layout
plt.savefig('UvonPhi.pdf')

#Fitten der Abstandsfunktion
def g(x, a, b):
    return b / x**a
def g1(x, m, b):
    return 1/m * 1/( x**2 ) + b
params1, covariance1 = curve_fit(g1, AbsCm, ULED)
errors1 = np.sqrt(np.diag(covariance1))
m = ufloat(params1[0], errors1[0])
b = ufloat(params1[1], errors1[1])
print(m)
print(b)
params2, covariance2 = curve_fit(g, AbsCm, ULED/Gain4)
errors2 = np.sqrt(np.diag(covariance2))
a = ufloat(params2[0], errors2[0])
b = ufloat(params2[1], errors2[1])
print(a)
print(b)
#
x1plot = np.linspace(10, 60)
plt.figure(1)
plt.title("Spannung als Funktion des Abstandes zwischen LED und PD")
plt.yscale("log")
plt.xscale("log")
plt.plot(AbsCm, ULED/Gain4, 'b+', label='Gefittete Messwerte')
plt.plot(AbsCmWeg, ULEDweg/Gain4, 'g+', label='Nicht gefittete Messwerte' )
#plt.plot(x1plot, g1(x1plot, *params1) , 'b-', label='Regression')
plt.plot(x1plot, g(x1plot, *params2) , 'g-', label='Regression 2')
plt.legend(loc="best")
plt.tight_layout
plt.savefig('Abst.pdf')

fig3, (ax1, ax2) = plt.subplots(2, 1)

ax1.set_title(r"Spannung als Funktion des Abstandes zwischen LED und PD")
ax1.set_xlabel("$r$/cm")
ax1.set_ylabel("$U(r)$/$\,$V")
ax1.plot(AbsCm, ULED/Gain4, 'b+', label='Gefittete Messwerte')
ax1.plot(AbsCmWeg, ULEDweg/Gain4, 'r+', label='Nicht gefittete Messwerte' )
ax1.plot(x1plot, g(x1plot, *params2) , 'g-', label='Regression')
ax1.legend(loc="best")

ax2.set_title(r"Spannung als Funktion des Abstandes zwischen LED und PD, doppelt logarithmisch aufgetragen")
ax2.set_xlabel("$r$/cm")
ax2.set_ylabel("$U(r)$/$\,$V")
ax2.plot(AbsCm, ULED/Gain4, 'b+', label='Gefittete Messwerte')
ax2.plot(AbsCmWeg, ULEDweg/Gain4, 'r+', label='Nicht gefittete Messwerte' )
ax2.plot(x1plot, g(x1plot, *params2) , 'g-', label='Regression')
ax2.set_yscale("log")
ax2.set_xscale("log")
ax2.legend(loc="best")

fig3.tight_layout()
plt.savefig('SubABst.pdf')
