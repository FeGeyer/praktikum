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
x20, y20 = np.genfromtxt('Steigung_T_20.txt', unpack=True)
x152, y152 = np.genfromtxt('Steigung_T_152.txt', unpack=True)
USkalaV, x20Skalamm, x152Skalamm = np.genfromtxt('xSkalenA.txt', unpack=True)
USkalaBV, xSkalaBmm = np.genfromtxt('xSkalaB.txt', unpack=True)
Abstaende = np.genfromtxt('Kurve.txt', unpack=True)
cx, cy = np.genfromtxt('ion.txt', unpack=True)
USkalaCV, xSkalaCmm = np.genfromtxt('xSkalaC.txt', unpack=True)

cy = cy + 5

USkalaBV = USkalaBV*5
USkalaCV = USkalaCV*5
#Spannungsmessgerät Kaputt, 5V sind nicht 5V, eher 6.5

def f (x, m, b):
    return m * x + b
paramsT20, covarianceT20 = curve_fit(f, x20Skalamm, USkalaV)
errorsT20 = np.sqrt(np.diag(covarianceT20))
mT20 = paramsT20[0]
bT20 = paramsT20[1]

paramsT152, covarianceT152 = curve_fit(f, x152Skalamm, USkalaV)
errorsT152 = np.sqrt(np.diag(covarianceT152))
mT152 = paramsT152[0]
bT152 = paramsT152[1]

print("Teil a)")
print("")
print("Umrechnungen Teil a):")
print("T=20")
print(mT20, "+/-" , errorsT20[0])
print(bT20, "+/-" , errorsT20[1])
print("T=152")
print(mT152, "+/-" , errorsT152[0])
print(bT152, "+/-" , errorsT152[1])


x20 = mT20*x20 + bT20
x152 = mT152*x152 + bT152

m20 = y20/x20
m152 = y152/x152

x20mm = np.array([0]*41)
x152mm = np.array([0]*41)

for i in range(41):
 x20mm[i] =i*5

for i in range(41):
 x152mm[i] = i*5

x20mm = mT20*x20mm + bT20
x152mm = mT152*x152mm + bT152

KDiff = 9-x20mm[31]

print("K aus Differntieller Verteilung: ", KDiff, "V")

plt.figure(1)
plt.xlim(0,11)
plt.ylim(-1, 120)
plt.xticks([0, 2, 4, 6, 8, 9, 10], ["0", "2", "4", "6", "8", r"$U_\mathrm{B}$", "10"])
plt.xlabel(r"$U/$V")
plt.ylabel(r"$f^\prime(U)/\mathrm{mm}\,\mathrm{V}^{-1}$")
plt.plot(x20mm, m20, 'rx', label=r"Diferentielle Energieverteilung $T=20 °C$")
plt.axvline(x=x20mm[31], ymin=0, ymax=1, color="k", ls=':', label=r"$U_{\mathrm{max}}$")
plt.axvline(x=9, ymin=0, ymax=1, color="k", ls='-.', label=r"$U_{\mathrm{B}}=9$V")
plt.annotate(s="K", xy=(x20mm[31]+(KDiff/2), 60), xytext=(9.5,70), arrowprops=dict(arrowstyle='->'))
plt.arrow(x20mm[31], 60, KDiff, 0)
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("T20.pdf")

plt.figure(2)
#plt.xlim(-0.015, 1)
#plt.ylim(0, 10.5)
plt.xlabel(r"$U/$V")
plt.ylabel(r"$f^\prime(U)/\mathrm{mm}\,\mathrm{V}^{-1}$")
plt.plot(x152mm, m152, 'rx', label=r"Diferentielle Energieverteilung $T=152 °C$")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("T152.pdf")

print("")
print("Teil b)")
print("")

paramsB, covarianceB = curve_fit(f, xSkalaBmm, USkalaBV)
errorsB = np.sqrt(np.diag(covarianceB))
mB = paramsB[0]
bB = paramsB[1]

print("Umrechnungen Teil b):")
print(mB, "+/-" , errorsB[0])
print(bB, "+/-" , errorsB[1])

Abstaende = mB * Abstaende + bB

KplusDU = mB * 22 + bB

DU = ufloat(np.mean(Abstaende), stats.sem(Abstaende))

KHertz = KplusDU - DU

print("E1-E0 = ",DU, "eV")
print("K2= ", KHertz, "V")

lamb = const.c/(DU * const.e / const.h)
print(lamb)

print("")
print("Teil c)")
print("")


paramsC, covarianceC = curve_fit(f, xSkalaCmm, USkalaCV)
errorsC = np.sqrt(np.diag(covarianceC))
mC = paramsC[0]
bC = paramsC[1]

print("Umrechnungen Teil c):")
print(mC, "+/-" , errorsC[0])
print(bC, "+/-" , errorsC[1])

paramsIon, covarianceIon = curve_fit(f, cx, cy)
errorsIon = np.sqrt(np.diag(covarianceIon))
mIon = ufloat(paramsIon[0], errorsIon[0])
bIon = ufloat(paramsIon[1], errorsIon[1])

EIon = -1*bIon/mIon
EIon = mC * EIon + bC

print("Ionisierungsenergie + K: " , EIon, "eV")
print("Ionisierungsenergie: " , EIon-(KDiff+KHertz)/2, "eV")

xSpace = np.linspace(90, 215)
plt.figure(3)
plt.xlim(unp.nominal_values(EIon)-0.25, 30.25)
plt.ylim(0, 180)
plt.xlabel(r"$U/$V")
plt.ylabel(r"$I(U)/\mathrm{mm}$")
plt.plot(mC*cx+bC, cy, 'rx', label=r"Ionisierungsmessung")
plt.plot(mC*xSpace+bC, f(xSpace, *paramsIon), 'r-', label=r"Regression")
plt.annotate(s=r"$U_{\mathrm{ion}} + \bar{K}$", xy=(unp.nominal_values(EIon), 0),
xytext=(16,10), arrowprops=dict(arrowstyle='-'))
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Ion.pdf")
