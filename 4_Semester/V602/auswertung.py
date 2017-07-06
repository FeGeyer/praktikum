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
d = 201.4 *10**(-12) #m
# Germanium
thetaGe1, IGe = np.genfromtxt('germanium.txt', unpack = True)
thetaGe1 = thetaGe1/2
EGe = const.h * const.c /(2*d*np.sin(2*np.pi*thetaGe1/360))/const.e
np.savetxt('germaniumE.txt', np.column_stack([EGe*10**(-3), IGe]),fmt="%.1f")

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

#Strontium
thetaSr1, ISr = np.genfromtxt('strontium.txt', unpack = True)
thetaSr1 = thetaSr1/2
ESr = const.h * const.c /(2*d*np.sin(2*np.pi*thetaSr1/360))/const.e
np.savetxt('strontiumE.txt', np.column_stack([ESr*10**(-3), ISr]),fmt="%.1f")

plt.figure(3)
ax = plt.gca()
ax.invert_xaxis()
plt.xlim(19.8, 13.6)
plt.ylim(30, 170)
plt.xlabel(r"$E / \mathrm{keV}$")
plt.ylabel(r"$I / \mathrm{N/s}$")
ax.plot(ESr*10**(-3), ISr, 'r-', label="Messwerte")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("strontium.pdf")

#Brom
thetaBr1, IBr = np.genfromtxt('brom.txt', unpack = True)
thetaBr1 = thetaBr1/2
EBr = const.h * const.c /(2*d*np.sin(2*np.pi*thetaBr1/360))/const.e
np.savetxt('bromE.txt', np.column_stack([EBr*10**(-3), IBr]),fmt="%.1f")

plt.figure(4)
ax = plt.gca()
ax.invert_xaxis()
plt.xlim(16.2, 11.8)
plt.ylim(14, 48)
plt.xlabel(r"$E / \mathrm{keV}$")
plt.ylabel(r"$I / \mathrm{N/s}$")
ax.plot(EBr*10**(-3), IBr, 'r-', label="Messwerte")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("brom.pdf")

# Zink
thetaZn1, IZn = np.genfromtxt('Zink.txt', unpack = True)
thetaZn1 = thetaZn1/2
EZn = const.h * const.c /(2*d*np.sin(2*np.pi*thetaZn1/360))/const.e
np.savetxt('zinkE.txt', np.column_stack([EZn*10**(-3), IZn]),fmt="%.1f")

plt.figure(5)
ax = plt.gca()
ax.invert_xaxis()
plt.xlim(10.6, 8.5)
plt.ylim(50, 1200)
plt.xlabel(r"$E / \mathrm{keV}$")
plt.ylabel(r"$I / \mathrm{N/s}$")
ax.plot(EZn*10**(-3), IZn, 'r-', label="Messwerte")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("zink.pdf")

# Zirkonium
thetaZr1, IZr = np.genfromtxt('strontium.txt', unpack = True)
thetaZr1 = thetaZr1/2
EZr = const.h * const.c /(2*d*np.sin(2*np.pi*thetaZr1/360))/const.e
np.savetxt('strontiumE.txt', np.column_stack([EZr*10**(-3), IZr]),fmt="%.1f")

plt.figure(6)
ax = plt.gca()
ax.invert_xaxis()
plt.xlim(19.7,13.6)
plt.ylim(33, 166)
plt.xlabel(r"$E / \mathrm{keV}$")
plt.ylabel(r"$I / \mathrm{N/s}$")
ax.plot(EZr*10**(-3), IZr, 'r-', label="Messwerte")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("strontium.pdf")

# Emission
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
