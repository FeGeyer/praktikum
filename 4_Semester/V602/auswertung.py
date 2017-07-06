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

thetaCu = thetaCu/2
d = 201.4 *10**(-12) #m
maxEnergie10 = const.h * const.c /(2*d*np.sin(2*np.pi*5/360))/const.e*10**(-3)
#maxEnergie104 = const.h * const.c /(2*d*np.sin(2*np.pi*5.2/360))/const.e
print("")
print("Maximale Energie bei Grenzwinkel 10°: ", maxEnergie10)

# K-Kanten
# Germanium K Kante
thetaGe = 15.95
kEnergieGe = const.h * const.c /(2*d*np.sin(2*np.pi*thetaGe/360))/const.e*10**(-3)
print("")
print("K-Kante Energie Germanium: ", kEnergieGe)

# Strontium
thetaSt = 10.95
kEnergieSt = const.h * const.c /(2*d*np.sin(2*np.pi*thetaSt/360))/const.e*10**(-3)
print("")
print("K-Kante Energie Strontium: ", kEnergieSt)

# Brom
thetaBr = 13
kEnergieBr = const.h * const.c /(2*d*np.sin(2*np.pi*thetaBr/360))/const.e*10**(-3)
print("")
print("K-Kante Energie Brom: ", kEnergieBr)

#Zink
thetaZn = 18.7
kEnergieZn = const.h * const.c /(2*d*np.sin(2*np.pi*thetaZn/360))/const.e*10**(-3)
print("")
print("K-Kante Energie Zink: ", kEnergieZn)

# Zirkonium
thetaZr = 9.85
kEnergieZr = const.h * const.c /(2*d*np.sin(2*np.pi*thetaZr/360))/const.e*10**(-3)
print("")
print("K-Kante Energie Zirkonium: ", kEnergieZr)
