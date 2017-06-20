import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
import scipy.constants as const
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

#-------------------------------------------------------------------------------
# Aufgabenteil a)
#-------------------------------------------------------------------------------

# Einlesen
UB1V, IB1nA = np.genfromtxt('blau11.txt', unpack=True)
UB2V, IB2nA = np.genfromtxt('blau2.txt', unpack=True)
UGeV, IGenA = np.genfromtxt('gelb.txt', unpack=True)
UGrV, IGrnA = np.genfromtxt('gruen.txt', unpack=True)
UUVV, IUVnA = np.genfromtxt('uv.txt', unpack=True)

IB1nA = IB1nA * 10**(-9)
IB2nA = IB2nA * 10**(-9)
IGenA = IGenA * 10**(-9)
IGrnA = IGrnA * 10**(-9)
IUVnA = IUVnA * 10**(-9)

UB1V = UB1V * (-1)
UB2V = UB2V * (-1)
UGeV = UGeV * (-1)
UGrV = UGrV * (-1)
UUVV = UUVV * (-1)

def f(U, m, n):
    return U*m + n

paramsB1, covarianceB1 = curve_fit(f, UB1V[2:9], np.sqrt(IB1nA[2:9]))
errorsB1 = np.sqrt(np.diag(covarianceB1))
mB1 = ufloat(paramsB1[0], errorsB1[0])
nB1 = ufloat(paramsB1[1], errorsB1[1])

paramsB2, covarianceB2 = curve_fit(f, UB2V[1:8], np.sqrt(IB2nA[1:8]))
errorsB2 = np.sqrt(np.diag(covarianceB2))
mB2 = ufloat(paramsB2[0], errorsB2[0])
nB2 = ufloat(paramsB2[1], errorsB2[1])

paramsGe, covarianceGe = curve_fit(f, UGeV[2:7], np.sqrt(IGenA[2:7]))
errorsGe = np.sqrt(np.diag(covarianceGe))
mGe = ufloat(paramsGe[0], errorsGe[0])
nGe = ufloat(paramsGe[1], errorsGe[1])

paramsGr, covarianceGr = curve_fit(f, UGrV[1:10], np.sqrt(IGrnA[1:10]))
errorsGr = np.sqrt(np.diag(covarianceGr))
mGr = ufloat(paramsGr[0], errorsGr[0])
nGr = ufloat(paramsGr[1], errorsGr[1])

paramsUV, covarianceUV = curve_fit(f, UUVV[2:7], np.sqrt(IUVnA[2:7]))
errorsUV = np.sqrt(np.diag(covarianceUV))
mUV = ufloat(paramsUV[0], errorsUV[0])
nUV = ufloat(paramsUV[1], errorsUV[1])

print("Ug:")
UgB1 = (-1) * nB1 / mB1
UgB2 = (-1) * nB2 / mB2
UgGe = (-1) * nGe / mGe
UgGr = (-1) * nGr / mGr
UgUV = (-1) * nUV / mUV
print("Gelb: ", UgGe , mGe, nGe)
print("Grün: ", UgGr , mGr, nGr)
print("Blau1: ", UgB2, mB2, nB2)
print("Blau2: ", UgB1, mB1, nB2)
print("UV: ", UgUV, mUV, nUV)

# Plot
x1plot = np.linspace(-0.2, 0.85)

plt.figure(1)
plt.yticks([0, 0.000005, 0.000010, 0.000015],["5", "10", "15"])
plt.ylim(0, 0.000015)
plt.xlim(-0.2, 0.85)
plt.xlabel(r"$U / \mathrm{V}$")
plt.ylabel(r"$\sqrt{I / \mathrm{nA}} \cdot 10^{6}$")
plt.plot(UGeV, np.sqrt(IGenA),'y.', label="")
plt.plot(UGrV, np.sqrt(IGrnA),'g.', label="")
plt.plot(UB2V, np.sqrt(IB2nA),'k.', label="")
plt.plot(UB1V, np.sqrt(IB1nA),'b.', label="")
plt.plot(UUVV, np.sqrt(IUVnA),'r.', label="")

plt.plot(x1plot, f(x1plot, *paramsGe), 'y--', label="Orange")
plt.plot(x1plot, f(x1plot, *paramsGr), 'g-.', label="Grün")
plt.plot(x1plot, f(x1plot, *paramsB1), 'b-', label="Violett 1")
plt.plot(x1plot, f(x1plot, *paramsB2), 'k--', label="Violett 2")
plt.plot(x1plot, f(x1plot, *paramsUV), 'r:', label="Violett 3")

#plt.plot(UB1V,np.sqrt(IB1nA),'b-.', label="")
#plt.plot(UB2V,np.sqrt(IB2nA),'k-.', label="")
#plt.plot(UGeV,np.sqrt(IGenA),'y-.', label="")
#plt.plot(UGrV,np.sqrt(IGrnA),'g-.', label="")
#plt.plot(UUVV,np.sqrt(IUVnA),'r-.', label="")

plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Photostrom.pdf")

#-------------------------------------------------------------------------------
# Aufgabenteil b)
#-------------------------------------------------------------------------------

#[B1, B2, Ge, Gr, UV]
lam = np.array([436, 434, 579, 546, 405])
lam = lam * 10 ** (-9)
nu = const.c/lam

print(nu)

Ug = unp.nominal_values(np.array([UgB1, UgB2, UgGe, UgGr, UgUV]))
Ugerr = unp.std_devs(np.array([UgB1, UgB2, UgGe, UgGr, UgUV]))

paramsAus, covarianceAus = curve_fit(f, nu, Ug, sigma=Ugerr)
errorsAus = np.sqrt(np.diag(covarianceAus))
mAus = ufloat(paramsAus[0], errorsAus[0])
nAus = ufloat(paramsAus[1], errorsAus[1])

print("h/e: ", "1. Theorie: ", const.h/const.e, "  |  2. Regression: ", mAus)
print("Austrittsarbeit: ", np.absolute(nAus), "eV")

#nu = nu/(10**15)

x2plot = np.linspace(0, 0.9)
x2plot = x2plot * 10**15
plt.figure(2)
plt.xlim(0*10**15, 0.9 * 10**15)
plt.ylabel(r"$U_g / \mathrm{eV}$")
plt.xlabel(r"$\nu / \mathrm{Hz}$")
plt.errorbar(nu, Ug, Ugerr,  fmt='r.', label="Wertepaare")
plt.plot(x2plot, f(x2plot, *paramsAus), 'r', label="")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Austritt.pdf")

#-------------------------------------------------------------------------------
# Aufgabenteil c)
#-------------------------------------------------------------------------------

Uc, IcnA = np.genfromtxt('kennline.txt', unpack=True)
plt.figure(3)
plt.xlabel(r"$U/ \mathrm{eV}$")
plt.ylabel(r"$I / \mathrm{nA}$")
plt.plot(Uc, IcnA, "bx", label="Messwerte")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Kennlinie.pdf")
