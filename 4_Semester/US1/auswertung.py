import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

# Impuls-Echo-Verfahren
# Einlesen der Daten
IlaengeMM, IampliV, ItMikros, IdeltatMikros, ITGCdB = np.genfromtxt('impuls.txt', unpack=True)
DlaengeMM, DampliV, DtMikros, DTGCdB = np.genfromtxt('durchschall.txt', unpack=True)
Apeak, AtMikros, AdeltatMikros = np.genfromtxt('auge.txt', unpack=True)
Cpeak, CtMikros = np.genfromtxt('cepstrum.txt', unpack=True)

IlaengeMM = unp.uarray(IlaengeMM, 0.02)
DlaengeMM = unp.uarray(DlaengeMM, 0.02)

Ic = (IlaengeMM*10**(-3)) / (0.5*IdeltatMikros*10**(-6))
Dc = (DlaengeMM*10**(-3)) / (DtMikros*10**(-6))

print(Ic)
print(np.mean(unp.nominal_values(Ic)))
print(np.sqrt(1/7*np.sum((unp.std_devs(Ic))**2)))

Icmean = ufloat(np.mean(unp.nominal_values(Ic)), stats.sem(unp.nominal_values(Ic)))
Dcmean = ufloat(np.mean(unp.nominal_values(Dc)), stats.sem(unp.nominal_values(Dc)))

#Lineare Regression
def f(x, m, n):
    return x/m + n
paramsI, covarianceI = curve_fit(f, unp.nominal_values(IlaengeMM)*10**(-3), 0.5*IdeltatMikros*10**(-6))
errorsI = np.sqrt(np.diag(covarianceI))
dtI = ufloat(paramsI[1], errorsI[1])
cI = ufloat(paramsI[0], errorsI[0])

print(cI, dtI)

paramsD, covarianceD = curve_fit(f, unp.nominal_values(DlaengeMM)*10**(-3), DtMikros*10**(-6))
errorsD = np.sqrt(np.diag(covarianceD))
dtD = ufloat(paramsD[1], errorsD[1])
cD = ufloat(paramsD[0], errorsD[0])

print(cD, dtD)

x1plot = np.linspace(0, 0.125)
plt.figure(1)
plt.grid()
plt.xlim(0, 0.125)
plt.xticks([0, 0.025, 0.05, 0.075, 0.1, 0.125],["0", "25", "50", "75", "100", "125"])
plt.ylim(0, 0.000045)
plt.yticks([0, 0.000005, 0.00001, 0.000015, 0.00002, 0.000025, 0.00003, 0.000035,
0.00004, 0.000045, 0.00005],["0", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50"])
plt.xlabel(r"$l/\mathrm{mm}$")
plt.ylabel(r"$t/\mathrm{\mu s}$")
plt.plot(unp.nominal_values(IlaengeMM*10**(-3)), 0.5*IdeltatMikros*10**(-6), 'bo', label='Messwerte')
plt.plot(x1plot, f(x1plot, *paramsI) , 'r-', label='Regression')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Imp.pdf")

plt.figure(2)
plt.grid()
plt.xlim(0, 0.125)
plt.xticks([0, 0.025, 0.05, 0.075, 0.1, 0.125],["0", "25", "50", "75", "100", "125"])
plt.ylim(0, 0.00005)
plt.yticks([0, 0.000005, 0.00001, 0.000015, 0.00002, 0.000025, 0.00003, 0.000035,
0.00004, 0.000045, 0.00005],["0", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50"])
plt.xlabel(r"$l/\mathrm{mm}$")
plt.ylabel(r"$t/\mathrm{\mu s}$")
plt.plot(unp.nominal_values(DlaengeMM*10**(-3)), DtMikros*10**(-6), 'bo', label='Messwerte')
plt.plot(x1plot, f(x1plot, *paramsD) , 'r-', label='Regression')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Dur.pdf")

IlaengeMM, IampliV, ItMikros, IdeltatMikros, ITGCdB = np.genfromtxt('ImpulsDump.txt', unpack=True)
DlaengeMM, DampliV, DtMikros, DTGCdB = np.genfromtxt('durchschallDump.txt', unpack=True)

DTGC = 10**(DTGCdB/20)
DGain = 10**(10/20)
DUber = DampliV/DTGC
DUber = DUber/DGain

ITGC = 10**(ITGCdB/20)
IGain = 10**(10/20)
IUber = IampliV/ITGC



def g(x, U, a):
    return U*np.exp(x * a)
paramsDD, covarianceDD = curve_fit(g, unp.nominal_values(DlaengeMM)*10**(-3), DUber)
errorsDD = np.sqrt(np.diag(covarianceDD))
DUnull = ufloat(paramsDD[0], errorsDD[0])
Dalpha = ufloat(paramsDD[1], errorsDD[1])

paramsDI, covarianceDI = curve_fit(g, 2*unp.nominal_values(IlaengeMM)*10**(-3), IUber)
errorsDI = np.sqrt(np.diag(covarianceDI))
IUnull = ufloat(paramsDI[0], errorsDI[0])
Ialpha = ufloat(paramsDI[1], errorsDI[1])

plt.figure(3)
plt.grid()
plt.xlim(0, 0.125)
plt.xticks([0, 0.025, 0.05, 0.075, 0.1, 0.125],["0", "25", "50", "75", "100", "125"])
plt.ylim(0, 1)
plt.xlabel(r"$l/\mathrm{mm}$")
plt.ylabel(r"$U/\mathrm{V}$")
plt.plot(unp.nominal_values(DlaengeMM*10**(-3)), DUber, 'bo', label='Messwerte')
plt.plot(x1plot, g(x1plot, *paramsDD) , 'r-', label='Regression')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("DurDump.pdf")

plt.figure(4)
x1plot = np.linspace(0, 0.250)
plt.grid()
plt.xlim(0, 0.250)
plt.xticks([0, 0.05, 0.1, 0.15, 0.2, 0.25 ],["0", "50", "100", "150", "200", "250"])
plt.ylim(0, 1)
plt.xlabel(r"$l/\mathrm{mm}$")
plt.ylabel(r"$U/\mathrm{V}$")
plt.plot(unp.nominal_values(2*IlaengeMM*10**(-3)), IUber, 'bo', label='Messwerte')
plt.plot(x1plot, g(x1plot, *paramsDI) , 'r-', label='Regression')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("ImpDump.pdf")

print(DUnull, Dalpha)
print(IUnull, Ialpha)

print(cD * (4.43*0.5) * 10 ** (-3))
print(cD * (8.67*0.5) * 10 ** (-3))
print(cD * (4.52*0.5) * 10 ** (-3))
print(cD * (8.76*0.5) * 10 ** (-3))
