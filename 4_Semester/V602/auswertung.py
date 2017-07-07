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
ax.axvline(x = 11.1670586743, ymin=0, ymax=1, ls='-.', color="k", label=r"$E_\mathrm{K}$")
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
ax.axvline(x = 16.2043769029, ymin=0, ymax=1, ls='-.', color="k", label=r"$E_\mathrm{K}$")
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
ax.axvline(x = 13.6317150382, ymin=0, ymax=1, ls='-.', color="k", label=r"$E_\mathrm{K}$")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("brom.pdf")

# Zink
thetaZn1, IZn = np.genfromtxt('Zink.txt', unpack = True)
thetaZn1 = thetaZn1/2
EZn = const.h * const.c /(2*d*np.sin(2*np.pi*thetaZn1/360))/const.e
np.savetxt('zinkE.txt', np.column_stack([EZn*10**(-3), IZn]),fmt="%.1f")

plt.figure(5)
plt.axes([0.1, 0.1, 0.8, 0.8])
ax = plt.gca()
ax.invert_xaxis()
plt.xlim(10.55, 8.6)
plt.ylim(50, 1200)
plt.xlabel(r"$E / \mathrm{keV}$")
plt.ylabel(r"$I / \mathrm{N/s}$")
ax.plot(EZn*10**(-3), IZn, 'r-', label="Messwerte")
ax.axvline(x = 9.60054213527, ymin=0, ymax=1, ls='-.', color="k", label=r"$E_\mathrm{K}$")
plt.legend(loc="best")
plt.axes([0.32975, 0.25, 0.3, 0.3])
ax = plt.gca()
ax.invert_xaxis()
plt.xlim(10, 9.25)
plt.xticks([9.25, 9.5, 9.75, 10], ["9.25", "9.50", "9.75", "10.00"])
plt.ylim(55, 110)
ax.plot(EZn[10:30] *10**(-3), IZn[10:30], 'r-')
ax.axvline(x = 9.60054213527, ymin=0, ymax=1, ls='-.', color="k")
plt.savefig("zink.pdf")

# Zirkonium
thetaZr1, IZr = np.genfromtxt('zirkonium.txt', unpack = True)
thetaZr1 = thetaZr1/2
EZr = const.h * const.c /(2*d*np.sin(2*np.pi*thetaZr1/360))/const.e
np.savetxt('zirkoniumE.txt', np.column_stack([EZr*10**(-3), IZr]),fmt="%.1f")

plt.figure(6)
ax = plt.gca()
ax.invert_xaxis()
plt.xlim(19.7,13.6)
plt.ylim(33, 166)
plt.xlabel(r"$E / \mathrm{keV}$")
plt.ylabel(r"$I / \mathrm{N/s}$")
ax.plot(EZr*10**(-3), IZr, 'r-', label="Messwerte")
ax.axvline(x = 17.9930435103, ymin=0, ymax=1, ls='-.', color="k", label=r"$E_\mathrm{K}$")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("zirkonium.pdf")

# Emission
thetaCu = thetaCu/2
maxEnergie10 = const.h * const.c /(2*d*np.sin(2*np.pi*5/360))/const.e
ECu = const.h * const.c /(2*d*np.sin(2*np.pi*thetaCu/360))/const.e
np.savetxt('emissionE.txt', np.column_stack([ECu*10**(-3), ICu]),fmt="%.3f")

plt.figure(1)
ax = plt.gca()
ax.invert_xaxis()
#plt.xlim(19.8, 11.8)
#plt.ylim(91, 141)
plt.xlabel(r"$E / \mathrm{keV}$")
plt.ylabel(r"$I / \mathrm{N/s}$")
ax.plot(ECu*10**(-3), ICu, 'r-', label="Messwerte")
#ax.axvline(x = 13.25, ymin=0, ymax=1, ls='-.', color="k", label=r"$\mathrm{K_{\alpha}}$")
#ax.axvline(x = 15.3, ymin=0, ymax=1, ls=':', color="k", label=r"$\mathrm{K_{\beta}}$")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("emission.pdf")
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

#Moseley
Ek = [kEnergieZn, kEnergieGe, kEnergieBr, kEnergieSr, kEnergieZr]
Z = [30, 32, 35, 38, 40]

#Fit
def f(x, m, b):
    return m * x + b

paramsR, covarianceR = curve_fit(f, Z, np.sqrt(Ek), p0 = [13.6, 0])
errorsR = np.sqrt(np.diag(covarianceR))
R = ufloat(paramsR[0], errorsR[0])
b = ufloat(paramsR[1], errorsR[1])

R = R**2
print("")
print("--------------------------------------")
print("Rydberg-Energie: ", R)

plt.figure(7)
x = np.linspace(29, 41)
plt.xlim(29, 41)
#plt.ylim(33, 166)
plt.xlabel(r"$Z$")
plt.ylabel(r"$E_k^{1/2} / \mathrm{eV}$")
plt.plot(Z, np.sqrt(Ek), 'r.', label="Messwerte")
plt.plot(x, f(x, *paramsR), 'b', label="Regression")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("moseley.pdf")

#Wismut
thetaWi, IWi = np.genfromtxt('wismut.txt', unpack = True)
thetaWi = thetaWi/2
EWi = const.h * const.c /(2*d*np.sin(2*np.pi*thetaWi/360))/const.e
np.savetxt('wismutE.txt', np.column_stack([EWi*10**(-3), IWi]),fmt="%.1f")

plt.figure(8)
ax = plt.gca()
ax.invert_xaxis()
plt.xlim(19.8, 11.8)
plt.ylim(91, 141)
plt.xlabel(r"$E / \mathrm{keV}$")
plt.ylabel(r"$I / \mathrm{N/s}$")
ax.plot(EWi*10**(-3), IWi, 'r-', label="Messwerte")
ax.axvline(x = 13.25, ymin=0, ymax=1, ls='-.', color="k", label=r"$\mathrm{K_{\alpha}}$")
ax.axvline(x = 15.3, ymin=0, ymax=1, ls=':', color="k", label=r"$\mathrm{K_{\beta}}$")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("wismut.pdf")

E1 = 13.25*10**(3)
E2 = 15.3*10**(3)
deltaE = E2- E1
sigmaWi = 83 - np.sqrt(4/const.alpha * np.sqrt(deltaE/13.6)- 5 * deltaE/13.6) * np.sqrt(1 + 19/32 * const.alpha**2 * deltaE/13.6)
print("")
print("Abschirmkonstante Wismut: ", sigmaWi)

deltaE2 = 15.71*10**3-13.42*10**(3)
sigmaWi2 = 83 - np.sqrt(4/const.alpha * np.sqrt(deltaE2/13.6)- 5 * deltaE2/13.6) * np.sqrt(1 + 19/32 * const.alpha**2 * deltaE2/13.6)
print("Literaturwert Wismut: ", sigmaWi2)

# Abschirmkonstanten für Emissionsspektrum
Kalpha = 8.077 *10**3
Kbeta = 9 *10**3
deltaK = Kbeta-Kalpha
sigmaK = 29 - np.sqrt(Kbeta/13.6)
sigmaL = 29 - np.sqrt(4*(deltaK)/13.6)
print("")
print("Sigma k: ", sigmaK)
print("Sigma L: ", sigmaL)

#Fit für HAlbwertsbreite
def f(x, m, b):
    return m * x + b

paramsKbeta, covarianceKbeta = curve_fit(f, ECu[78:80], ICu[78:80])
errorsKbeta = np.sqrt(np.diag(covarianceKbeta))
mKbeta = ufloat(paramsR[0], errorsR[0])
bKbeta = ufloat(paramsR[1], errorsR[1])
#print(mKbeta)
#print(bKbeta)
