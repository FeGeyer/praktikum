import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

# Erster Teil
# Einlesen der Daten
pumpPro15, fmaxHz15, fmeanHz15, inten15 = np.genfromtxt('winkel15.txt', unpack=True)
pumpPro30, fmaxHz30, fmeanHz30, inten30 = np.genfromtxt('winkel30.txt', unpack=True)
pumpPro60, fmaxHz60, fmeanHz60, inten60 = np.genfromtxt('winkel60.txt', unpack=True)

# Bestimmung der Momentangeschwindigkeit
f0Hz = 2 * 10 ** 6
c = 1800
v15 = fmeanHz15*c/(2*f0Hz) * 1/np.cos(80.06)
v30 = fmeanHz30*c/(2*f0Hz) * 1/np.cos(70.53)
v60 = fmeanHz60*c/(2*f0Hz) * 1/np.cos(54.74)
print(1/np.cos(80.06))
print(pumpPro60)
print(fmeanHz15)
print(v15)
print(pumpPro60)
print(fmeanHz30)
print(v30)
print(pumpPro60)
print(fmeanHz60)
print(v60)

# Plots
plt.figure(1)
plt.xlabel(r"$v/ \mathrm{\frac{m}{s}}$")
plt.ylabel(r"$\frac{\mathrm{\Delta} \nu}{\mathrm{cos} \alpha} / \mathrm{kHz}$")
plt.plot(v15, fmeanHz15/np.cos(80.06) * 10 ** (-3), 'ro', label="Dopplerwinkel 80.06°")
plt.ylim(0.5, 8)
plt.xlim(0.5, 3.5)
plt.legend(loc="best")
plt.tight_layout
plt.savefig('a15.pdf')

plt.figure(2)
plt.xlabel(r"$v/ \mathrm{\frac{m}{s}}$")
plt.ylabel(r"$\frac{\mathrm{\Delta} \nu}{\mathrm{cos} \alpha} / \mathrm{kHz}$")
plt.plot(v30, fmeanHz30/np.cos(70.53) * 10 ** (-3), 'go', label="Dopplerwinkel 70.53°")
plt.ylim(0.4, 4)
plt.xlim(0.4, 1.8)
#plt.xticks([0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6],
#            [])
plt.legend(loc="best")
plt.tight_layout
plt.savefig('a30.pdf')

plt.figure(3)
plt.xlabel(r"$v/ \mathrm{\frac{m}{s}}$")
plt.ylabel(r"$\frac{\mathrm{\Delta} \nu}{\mathrm{cos} \alpha} / \mathrm{kHz}$")
plt.plot(v60, fmeanHz60/np.cos(54.74) * 10 ** (-3), 'bo', label="Dopplerwinkel 54.74°")
plt.ylim(0.4, 4.5)
plt.xlim(0.4, 2)
plt.legend(loc="best")
plt.tight_layout
plt.savefig('a60.pdf')

# Aufgabenteil 2
# Auslesen der Daten
tiefeMikroS45, fmaxHz45, fmeanHz45, inten45 = np.genfromtxt('aufgabe2_45.txt', unpack=True)
tiefeMikroS70, fmaxHz70, fmeanHz70, inten70 = np.genfromtxt('aufgabe2_70.txt', unpack=True)

tiefeMM45 = tiefeMikroS45 * 1.5
tiefeMM70 = tiefeMikroS70 * 1.5

# Bestimmung der Momentangeschwindigkeit
v45 = fmeanHz45*c/(2*f0Hz) * 1/np.cos(80.06)
v70 = fmeanHz70*c/(2*f0Hz) * 1/np.cos(80.06)

# Plots
plt.figure(4)
plt.xlabel(r"$\mathrm{Messtiefe} / \mathrm{mm}$")
plt.ylabel(r"$\mathrm{Streuintensität} / \mathrm{\frac{V^2}{s} * 10}$")
plt.plot(tiefeMM45, inten45 *10 ** (-2), 'ro', label="Streuintensität")
plt.plot(tiefeMM45, v45, 'bo', label="Momentangeschwindigkeit")
plt.xlim(12, 25)
plt.legend(loc="best")
plt.tight_layout
plt.savefig('b45.pdf')

plt.figure(5)
plt.xlabel(r"$\mathrm{Messtiefe} / \mathrm{mm}$")
plt.ylabel(r"$\mathrm{Streuintensität} / \mathrm{\frac{V^2}{s} * 10}$")
plt.plot(tiefeMM70, inten70 *10 ** (-2), 'ro', label="Streuintensität")
plt.plot(tiefeMM70, v70, 'bo', label="Momentangeschwindigkeit")
plt.xlim(12, 25)
plt.ylim(0, 8)
plt.legend(loc="best")
plt.tight_layout
plt.savefig('b70.pdf')
