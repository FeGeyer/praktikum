import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
import scipy.constants as const
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

# Emissionsspektrum
# Einlesen der Daten
thetaCu, ICu = np.genfromtxt('emission.txt', unpack = True)
plt.figure(1)
#plt.xlim(29.5, 40.5)
plt.xlabel(r"$2 \Theta$")
plt.ylabel(r"$I / \mathrm{N/s}$")
plt.plot(thetaCu, ICu, 'r-', label="Messwerte")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("emission.pdf")

# Einlesen und in Energie umwandeln
# Germanium
d = 201.4 *10**(-12) #m
thetaGe1, IGe = np.genfromtxt('germanium.txt', unpack = True)
thetaGe1 = thetaGe1/2
EGe = const.h * const.c /(2*d*np.sin(2*np.pi*thetaGe1/360))/const.e


plt.figure(2)
ax = plt.gca()
ax.invert_xaxis()
plt.xlim(12.5, 9.9)
plt.ylim(30.5, 70.5)
plt.xlabel(r"$E / \mathrm{keV}$")
plt.ylabel(r"$I / \mathrm{N/s}$")
ax.plot(EGe*10**(-3), IGe, 'r-', label="Messwerte")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("germanium.pdf")

np.savetxt('germaniumE.txt', np.column_stack([EGe*10**(-3), IGe]),fmt="%.3f")

thetaCu = thetaCu/2
maxEnergie10 = const.h * const.c /(2*d*np.sin(2*np.pi*5/360))/const.e
#maxEnergie104 = const.h * const.c /(2*d*np.sin(2*np.pi*5.2/360))/const.e
print("")
print("Maximale Energie bei Grenzwinkel 10°: ", maxEnergie10*10**(-3))

# K-Kanten
# Germanium K Kante
thetaGe = 16
kEnergieGe = const.h * const.c /(2*d*np.sin(2*np.pi*thetaGe/360))/const.e
print("")
print("K-Kante Energie Germanium: ", kEnergieGe*10**(-3))

# Strontium
thetaSt = 10.95
kEnergieSr = const.h * const.c /(2*d*np.sin(2*np.pi*thetaSt/360))/const.e
print("")
print("K-Kante Energie Strontium: ", kEnergieSr*10**(-3))

# Brom
thetaBr = 13.05
kEnergieBr = const.h * const.c /(2*d*np.sin(2*np.pi*thetaBr/360))/const.e
print("")
print("K-Kante Energie Brom: ", kEnergieBr*10**(-3))

#Zink
thetaZn = 18.7
kEnergieZn = const.h * const.c /(2*d*np.sin(2*np.pi*thetaZn/360))/const.e
print("")
print("K-Kante Energie Zink: ", kEnergieZn*10**(-3))

# Zirkonium
thetaZr = 9.85
kEnergieZr = const.h * const.c /(2*d*np.sin(2*np.pi*thetaZr/360))/const.e
print("")
print("K-Kante Energie Zirkonium: ", kEnergieZr*10**(-3))

# Abschirmkonstanten
# Germanium
z = 32
sigmaGe = z - np.sqrt(kEnergieGe/13.6-(const.alpha**2 * z**4)/4)
print("")
print("-------------------------------------------")
print("")
print("Abschirmkonstante Germanium: ", sigmaGe)

#Strontium
z = 38
sigmaSr = z - np.sqrt(kEnergieSr/13.6-(const.alpha**2 * z**4)/4)
print("")
print("Abschirmkonstante Strontium: ", sigmaSr)

# Brom
z = 35
sigmaBr = z - np.sqrt(kEnergieBr/13.6-(const.alpha**2 * z**4)/4)
print("")
print("Abschirmkonstante Brom: ", sigmaBr)

# Zink
z = 30
sigmaZn = z - np.sqrt(kEnergieZn/13.6-(const.alpha**2 * z**4)/4)
print("")
print("Abschirmkonstante Zink: ", sigmaZn)

# Zirkonium
z = 40
sigmaZr = z - np.sqrt(kEnergieZr/13.6-(const.alpha**2 * z**4)/4)
print("")
print("Abschirmkonstante Zirkonium: ", sigmaZr)
