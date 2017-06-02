import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
import scipy.constants as const
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

# Einlesen
Ukenn, I1, I2, I3, I4, I5 = np.genfromtxt('kennlinien.txt', unpack=True)

# Plot
plt.figure(1)
plt.xlim(0, 260)
plt.xlabel(r"$U / \mathrm{V}$")
plt.ylabel(r"$I / \mathrm{mA}$")
plt.plot(Ukenn, I1,'.', label="Kennlinie 1")
plt.plot(Ukenn, I2,'.', label="Kennlinie 2")
plt.plot(Ukenn, I3,'.', label="Kennlinie 3")
plt.plot(Ukenn, I4,'.', label="Kennlinie 4")
plt.plot(Ukenn, I5,'.', label="Kennlinie 5")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("kennlinien.pdf")

# Sättigungsstrom
Is1 = np.max(I1) * 10**(-3)
Is2 = np.max(I2) * 10**(-3)
Is3 = np.max(I3) * 10**(-3)
Is4 = np.max(I4) * 10**(-3)
Is5 = np.max(I5) * 10**(-3)


print("")
print("Aufgabenteil a): ")
print(Is1)
print(Is2)
print(Is3)
print(Is4)
print(Is5)

# Aufgabenteil b)
# Werte logarithmieren
logU = np.log(Ukenn[1:7])
logI = np.log(I5[1:7])
def f (x, m, b):
    return m * x + b
paramsRaum, covarianceRaum = curve_fit(f, logU, logI, p0=[1.5, 2])
errorsRaum = np.sqrt(np.diag(covarianceRaum))
exponent = ufloat(paramsRaum[0], errorsRaum[0])
print("")
print("Aufgabenteil b): ")
print(exponent)


plt.figure(2)
x1 = np.linspace(2, logU[5])
plt.xlim(2, 5.5)
plt.xlabel(r"$\mathrm{ln}(U)$")
plt.ylabel(r"$\mathrm{ln}(I_5)$")
plt.plot(logU, logI,'r.', label="Gefittete Werte")
plt.plot(x1, f(x1, *paramsRaum), 'b-', label="Regression")
plt.axvline(x = logU[5], ymin=0, ymax=1, ls=':', color='k')
plt.plot(np.log(Ukenn[7:25]), np.log(I5[7:25]), 'g.', label="Nicht-gefittete Werte")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("raumladung.pdf")

# Aufgabenteil c)
U, I6 = np.genfromtxt('anlaufstrom.txt', unpack=True)
I6 = I6 * 10**(-9) # Umrechnen in Ampere
U = U - 1*10**6 * I6
logI6 = np.log(I6)
U1 = U[0:9]
logI61 = logI6[0:9]

paramsA, covarianceA = curve_fit(f, U1, logI61)
errorsA= np.sqrt(np.diag(covarianceA))
exponent2 = ufloat(paramsA[0], errorsA[0])

T = - (const.e)/(const.k*exponent2)
print("")
print("Aufgabenteil c): ")
print(T)

plt.figure(3)
x2 = np.linspace(-0.015, U1[8])
plt.xlim(-0.015, 1)
plt.ylim(-23, -18)
plt.xlabel(r"$U / \mathrm{V}$")
plt.ylabel(r"$\mathrm{ln}(I_6)$")
plt.plot(U1, logI61, 'r.', label="Gefittete Werte")
plt.axvline(x = U1[8], ymin=0, ymax=1, ls=':', color='k')
plt.plot(U[9:12], logI6[9:12], 'g.', label="Nicht-gefittete Werte")
plt.plot(x2, f(x2, *paramsA), 'b-', label="Regression")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("anlauf.pdf")

# Aufgabenteil d)
Uh, Ih = np.genfromtxt('leistung.txt', unpack=True)
NWL = 0.9
sigma = 5.7*10**(-12)
g = 0.32
eta = 0.28
T2 = ((Uh * Ih - NWL)/(sigma * g * eta))**(1/4)
T2mean = ufloat(np.mean(T2), stats.sem(T2))
print("")
print("Aufgabenteil d): ")
print(T2)
print(T2mean)

# Aufgabenteil e)
e0Phi1 = -np.log((Is1 * (const.h)**3)/(4 * g *10**(-4) * np.pi * const.e * const.m_e * (const.k)**2 * T2[0]**2)) * const.k * T2[0] /const.e
e0Phi2 = -np.log((Is2 * (const.h)**3)/(4 * g *10**(-4) * np.pi * const.e * const.m_e * (const.k)**2 * T2[1]**2)) * const.k * T2[1] /const.e
e0Phi3 = -np.log((Is3 * (const.h)**3)/(4 * g *10**(-4) * np.pi * const.e * const.m_e * (const.k)**2 * T2[2]**2)) * const.k * T2[2] /const.e
e0Phi4 = -np.log((Is4 * (const.h)**3)/(4 * g *10**(-4) * np.pi * const.e * const.m_e * (const.k)**2 * T2[3]**2)) * const.k * T2[3] /const.e
e0Phi5 = -np.log((Is5 * (const.h)**3)/(4 * g *10**(-4) * np.pi * const.e * const.m_e * (const.k)**2 * T2[4]**2)) * const.k * T2[4] /const.e

e0Phi = [e0Phi1, e0Phi2, e0Phi3, e0Phi4, e0Phi5]
e0Phimean = ufloat(np.mean(e0Phi), stats.sem(e0Phi))
print("")
print("Aufgabenteil e): ")
print(e0Phi)
print(e0Phimean)
