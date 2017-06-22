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
dyU = (dyU2 - dyU1)*10**(-3) # V

gdR = gdR1 - gdR2 # Ohm
gdU = (gdU2 - gdU1)*10**(-3) # V

ndR = ndR1 - ndR2 # Ohm
ndU = (ndU2 - ndU1)*10**(-3) # V

# Qreal ausrechnen
dyQ = dyM/(dyL * dyRho) * 100 # mm^2
gdQ = gdM/(gdL * gdRho) * 100 # mm^2
ndQ = ndM/(ndL * ndRho) * 100 # mm^2

# Mit Widerstandsmethode
dyChiR = 2 * dyR/R3 * F/(dyQ)
gdChiR = 2 * gdR/R3 * F/(gdQ)
ndChiR = 2 * ndR/R3 * F/(ndQ)

dyChiRmean = ufloat(np.mean(dyChiR), stats.sem(dyChiR))
gdChiRmean = ufloat(np.mean(gdChiR), stats.sem(gdChiR))
ndChiRmean = ufloat(np.mean(ndChiR), stats.sem(ndChiR))


print("")
print("Q von Nd: ", ndQ)
print("Q von Gd: ", gdQ)
print("Q von Dy: ", dyQ)
print("")
print("Chi von Dy mit R: ", dyChiRmean)
print("Chi von Gd mit R: ", gdChiRmean)
print("Chi von Nd mit R: ", ndChiRmean)

# Mit Spannungsmethode
dyChiU = 4 * F/dyQ * dyU/0.9
gdChiU = 4 * F/gdQ * gdU/0.9
ndChiU = 4 * F/ndQ * ndU/0.9

dyChiUmean = ufloat(np.mean(dyChiU), stats.sem(dyChiU))
gdChiUmean = ufloat(np.mean(gdChiU), stats.sem(gdChiU))
ndChiUmean = ufloat(np.mean(ndChiU), stats.sem(ndChiU))

dyUmean = ufloat(np.mean(dyU), stats.sem(dyU))
gdUmean = ufloat(np.mean(gdU), stats.sem(gdU))
ndUmean = ufloat(np.mean(ndU), stats.sem(ndU))

print("")
print("Delta U von Nd: ", ndUmean)
print("Delta U von Gd: ", gdUmean)
print("Delta U von Dy: ", dyUmean)
print("")
print("Chi von Dy mit U: ", dyChiUmean)
print("Chi von Gd mit U: ", gdChiUmean)
print("Chi von Nd mit U: ", ndChiUmean)

# Theoretische Bestimmung

# Konstanten
T = 293.15 # Kelvin
muB = 1/2 * const.e/const.m_e * const.hbar

ndMol = 336.4822 # g/mol
gdMol = 362.4982 # g/mol
dyMol = 372.9982 # g/mol

ndGj = 8/11
gdGj = 2
dyGj = 4/3

ndJ = 4.5
gdJ = 3.5
dyJ = 7.5

# Bestimmung von N

ndN = 2 * const.N_A * ndRho/ndMol * 100**3 # 1/m^3
gdN = 2 * const.N_A * gdRho/gdMol * 100**3 # 1/m^3
dyN = 2 * const.N_A * dyRho/dyMol * 100**3 # 1/m^3

# Bestimmung von Chi

ndChi = (const.mu_0 * muB**2 * ndGj**2 * ndN * ndJ*(ndJ+1))/(3 * const.k * T)
gdChi = (const.mu_0 * muB**2 * gdGj**2 * gdN * gdJ*(gdJ+1))/(3 * const.k * T)
dyChi = (const.mu_0 * muB**2 * dyGj**2 * dyN * dyJ*(dyJ+1))/(3 * const.k * T)

#ndChi = np.round(ndChi, 4)

#print("N für Nd: ", ndN)
#print("N für Gd: ", gdN)
#print("N für Dy: ", dyN)

dyRmean = ufloat(np.mean(dyR), stats.sem(dyR))
gdRmean = ufloat(np.mean(gdR), stats.sem(gdR))
ndRmean = ufloat(np.mean(ndR), stats.sem(ndR))

print("")
print("Delta R für Nd: ", ndRmean)
print("Delta R für Gd: ", gdRmean)
print("Delta R für Dy: ", dyRmean)

print("")
print("Theoretisches Chi von Dy: ", np.round(dyChi, 5))
print("Theoretisches Chi von Gd: ", np.round(gdChi, 4))
print("Theoretisches Chi von Nd: ", np.round(ndChi, 6))
