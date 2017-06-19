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

UB1V = UB1V * (-1)
UB2V = UB2V * (-1)
UGeV = UGeV * (-1)
UGrV = UGrV * (-1)
UUVV = UUVV * (-1)

def f(U, m, n):
    return U*m + n

paramsB1, covarianceB1 = curve_fit(f, np.sqrt(IB1nA), UB1V)
errorsB1 = np.sqrt(np.diag(covarianceB1))
mB1 = ufloat(paramsB1[0], errorsB1[0])
nB1 = ufloat(paramsB1[1], errorsB1[1])

paramsB2, covarianceB2 = curve_fit(f, np.sqrt(IB2nA), UB2V)
errorsB2 = np.sqrt(np.diag(covarianceB2))
mB2 = ufloat(paramsB2[0], errorsB2[0])
nB2 = ufloat(paramsB2[1], errorsB2[1])

paramsGe, covarianceGe = curve_fit(f, np.sqrt(IGenA), UGeV)
errorsGe = np.sqrt(np.diag(covarianceGe))
mGe = ufloat(paramsGe[0], errorsGe[0])
nGe = ufloat(paramsGe[1], errorsGe[1])

paramsGr, covarianceGr = curve_fit(f, np.sqrt(IGrnA), UGrV)
errorsGr = np.sqrt(np.diag(covarianceGr))
mGr = ufloat(paramsGr[0], errorsGr[0])
nGr = ufloat(paramsGr[1], errorsGr[1])

paramsUV, covarianceUV = curve_fit(f, np.sqrt(IUVnA), UUVV)
errorsUV = np.sqrt(np.diag(covarianceUV))
mUV = ufloat(paramsUV[0], errorsUV[0])
nUV = ufloat(paramsUV[1], errorsUV[1])

print("Ug:")
UgB1 = nB1
UgB2 = nB2
UgGe = nGe
UgGr = nGr
UgUV = nUV
print("Gelb: ", UgGe)
print("Grün: ", UgGr)
print("Blau1: ", UgB1)
print("Blau2: ", UgB2)
print("UV: ", UgUV)


# Plot
x1plot = np.linspace(-0.2, 0.85)

plt.figure(1)
plt.ylim(0, 0.75)
plt.xlim(-0.2, 0.4)
plt.xlabel(r"$U / \mathrm{V}$")
plt.ylabel(r"$I / \mathrm{mA}$")
plt.plot(np.sqrt(IB1nA),UB1V,'b.', label="Blau 1")
plt.plot(np.sqrt(IB2nA),UB2V,'k.', label="Blau 2")
plt.plot(np.sqrt(IGenA),UGeV,'y.', label="Gelb")
plt.plot(np.sqrt(IGrnA),UGrV,'g.', label="Grün")
plt.plot(np.sqrt(IUVnA),UUVV,'r.', label="UV")

plt.plot(x1plot, f(x1plot, *paramsB1), 'b', label="")
plt.plot(x1plot, f(x1plot, *paramsB2), 'k', label="")
plt.plot(x1plot, f(x1plot, *paramsGe), 'y', label="")
plt.plot(x1plot, f(x1plot, *paramsGr), 'g', label="")
plt.plot(x1plot, f(x1plot, *paramsUV), 'r', label="")

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
lam = np.array([407, 435, 578, 546, 365])
lam = lam * 10 ** (-9)
nu = const.c/lam

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
plt.errorbar(nu, Ug, Ugerr,  fmt='r.', label="Paare")
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
