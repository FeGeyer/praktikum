import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
import scipy.constants as const
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

# Teil a)
# Einlesen
f, Ukenn = np.genfromtxt('filterkurve.txt', unpack=True)

# Plot
plt.figure(1)
plt.xlim(29.5, 40.5)
plt.xlabel(r"$f / \mathrm{kHz}$")
plt.ylabel(r"$U / \mathrm{V}$")
plt.plot(f, Ukenn*10**(-3), 'r.', label="Filterkurve")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("filterkurve.pdf")

# Aufgabenteil b)
# Einlesen
dyU1, dyR1, dyR2, dyU2 = np.genfromtxt('dy.txt', unpack=True)
gdU1, gdR1, gdR2, gdU2 = np.genfromtxt('gd.txt', unpack=True)
ndU1, ndR1, ndR2, ndU2 = np.genfromtxt('nd.txt', unpack=True)

# Umrechnen von R
dyR1 = dyR1 * 5 * 10**(-3) # Ohm
gdR1 = gdR1 * 5 * 10**(-3) # Ohm
ndR1 = ndR1 * 5 * 10**(-3) # Ohm

dyR2 = dyR2 * 5 * 10**(-3) # Ohm
gdR2 = gdR2 * 5 * 10**(-3) # Ohm
ndR2 = ndR2 * 5 * 10**(-3) # Ohm

# Konstanten Materialien
dyRho = 7.8 #g/cm^3
gdRho = 7.40 #g/cm^3
ndRho = 7.24 #g/cm^3

dyL = 15.5 # cm
gdL = 16.5 # cm
ndL = 16 # cm

dyM = 15.1 # g
gdM = 14.08  # g
ndM = 9.0  # g

R3 = 998 # Ohm

# Konstanten Spule
n = 250
F = 86.6 # mm²
l = 135 # mm
R = 0.7 # ohm

# Differenzen bilden
dyR = dyR1 - dyR2 # Ohm
dyU = dyU2 - dyU1 # Ohm

gdR = gdR1 - gdR2 # Ohm
gdU = gdU2 - gdU1 # Ohm

ndR = ndR1 - ndR2 # Ohm
ndU = ndU2 - ndU1 # Ohm

# Qreal ausrechnen
dyQ = dyM/(dyL * dyRho) * 100 # mm^2
gdQ = gdM/(gdL * gdRho) * 100 # mm^2
ndQ = ndM/(ndL * ndRho) * 100 # mm^2

# Mit zweiter Methode
dyChi = 2 * dyR/R3 * F/(dyQ)
gdChi = 2 * gdR/R3 * F/(gdQ)
ndChi = 2 * ndR/R3 * F/(ndQ)

dyChimean = ufloat(np.mean(dyChi), stats.sem(dyChi))
gdChimean = ufloat(np.mean(gdChi), stats.sem(gdChi))
ndChimean = ufloat(np.mean(ndChi), stats.sem(ndChi))
print("Chi von Dy: ", dyChimean)
print("Chi von Gd: ", gdChimean)
print("Chi von Nd: ", ndChimean)
