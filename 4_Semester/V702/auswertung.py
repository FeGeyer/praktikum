import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
import scipy.constants as const
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

tInmin, NIn = np.genfromtxt('indium.txt', unpack=True)
nrraus, NAgraus = np.genfromtxt('silberaus.txt', unpack=True)
nr, NAg = np.genfromtxt('silber.txt', unpack=True)

NIn = NIn/240
NAg = NAg/8

N0 = 164/900

tIns = tInmin * 60
tAgs = nr*8

#Lineare Regression für Kupfer
def f(t, m, n):
    return t*m + n
paramsIn, covarianceIn = curve_fit(f, tIns, np.log(NIn-N0), sigma=np.sqrt(NIn)-np.sqrt(N0))
errorsIn = np.sqrt(np.diag(covarianceIn))
InMu = ufloat(paramsIn[0], errorsIn[0])
InN0 = ufloat(paramsIn[1], errorsIn[1])
InN0 = unp.exp(InN0)
print("Indium:")
print("mu:",InMu)
print("N0/s:",InN0)
print("T1/2 /s:", np.log(2)/-InMu)

x1plot = np.linspace(0, 62)
plt.figure(1)
plt.xlabel(r"$t/ \mathrm{min}$")
plt.ylabel(r"$N/ \mathrm{s}^{-1}$")
plt.yscale("log")
plt.xlim(0,62)
#plt.ylim(5,15)
#plt.plot(tInmin, NIn, "r.", label="Messwerte Indium")
plt.errorbar(tInmin, NIn-N0, np.sqrt(NIn)-np.sqrt(N0),  fmt='r.', label="Messwerte Indium")
plt.plot(x1plot, np.exp(f(x1plot*60, *paramsIn)), 'b', label="Regression")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Indium.pdf")

paramsAg1, covarianceAg1 = curve_fit(f, tAgs[0:5], (NAg[0:5]-N0), sigma=np.sqrt(NAg[0:5])-np.sqrt(N0))
errorsAg1 = np.sqrt(np.diag(covarianceAg1))
Ag1Mu = ufloat(paramsAg1[0], errorsAg1[0])
Ag1N0 = ufloat(paramsAg1[1], errorsAg1[1])
#Ag1N0 = unp.exp(Ag1N0)
print("Silber-kurz:")
print("mu:",Ag1Mu)
print("N0/s:",Ag1N0)
print("T1/2 /s:",(Ag1N0/2) * (-1/Ag1Mu))
print(np.log(2)/-Ag1Mu)

paramsAg2, covarianceAg2 = curve_fit(f, tAgs[12:35], (NAg[12:35]-N0), sigma=np.sqrt(NAg[12:35])-np.sqrt(N0))
errorsAg2 = np.sqrt(np.diag(covarianceAg2))
Ag2Mu = ufloat(paramsAg2[0], errorsAg2[0])
Ag2N0 = ufloat(paramsAg2[1], errorsAg2[1])
#Ag2N0 = unp.exp(Ag2N0)
print("Silber-lang:")
print("mu:",Ag2Mu)
print("N0/s:",Ag2N0)
print("T1/2 /s:",(Ag2N0/2) * (-1/Ag2Mu))
print(np.log(2)/-Ag2Mu)

x2plot = np.linspace(0,7)
plt.figure(2)
plt.xlabel(r"$t/ \mathrm{min}$")
plt.ylabel(r"$N/ \mathrm{s}^{-1}$")
#plt.xlim(0,8.2)
plt.ylim(0,23)
#plt.plot(tAgs/60, NAg, "r.", label="Messwerte Silber")
plt.errorbar(tAgs/60, NAg-N0, np.sqrt(NAg)-np.sqrt(N0),  fmt='r.', label="Messwerte Silber")
plt.plot(nrraus*8/60, NAgraus/8-N0, 'kx', label="Nicht beachtete Messwerte")
plt.plot(x2plot, (f(x2plot*60, *paramsAg1)), 'y', label="Regression")
plt.plot(x2plot, (f(x2plot*60, *paramsAg2)), 'b', label="Regression")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Silber.pdf")

#np.savetxt('SilberTabelle.txt', np.column_stack([nr[0:14]*8, NAg[0:14]*8, nr[14:28]*8,
#NAg[14:28]*8, nr[28:42]*8, NAg[28:42]*8, nr[42:56]*8, NAg[42:56]*8 ]), fmt="%.0f" )
