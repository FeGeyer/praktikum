import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
import scipy.constants as const
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

# 7a
# Einlesen der Daten
# Kupfer
CuD, CuT, CuN = np.genfromtxt('CUGamma.txt', unpack=True)
# Nullrate abziehen
Nulleffekt = 1041/1000
# In SI umrechnen
CuD = CuD *1*10**(-3)
#Cu = unp.uarray(CuN, np.sqrt(CuN)/CuT)
# logarithmiert
ln_Cu = np.log(CuN/CuT - Nulleffekt)

# Zink
ZnD, ZnT, ZnN = np.genfromtxt('ZNGamma.txt', unpack=True)
ZnD = ZnD *1*10**(-3)
Zn = unp.uarray(ZnN, np.sqrt(ZnN)/ZnT)
ln_Zn = np.log(ZnN/ZnT - Nulleffekt)

#Lineare Regression für Kupfer
def f(x, m, n):
    return x*(-m) + n
paramsCu, covarianceCu = curve_fit(f, CuD, ln_Cu)
errorsCu = np.sqrt(np.diag(covarianceCu))
CuN0 = ufloat(paramsCu[1], errorsCu[1])
CuN0 = unp.exp(CuN0)
CuMu = ufloat(paramsCu[0], errorsCu[0])


#Lineare Regression für Zink
paramsZn, covarianceZn = curve_fit(f, ZnD, ln_Zn)
errorsZn = np.sqrt(np.diag(covarianceZn))
ZnN0 = ufloat(paramsZn[1], errorsZn[1])
ZnN0 = unp.exp(ZnN0)
ZnMu = ufloat(paramsZn[0], errorsZn[0])

print("")
print("Aufgabenteil a)")
print("Kupfer (N(0), mu): ", CuN0, CuMu)
print("Zink (N(0), mu): ",ZnN0, ZnMu)

# Plots
# Kupfer

x1plot = np.linspace(0, 18*10**(-3))
plt.figure(1)
plt.xlabel(r"$D/ \mathrm{mm}$")
plt.ylabel(r"$N / \mathrm{s}^{-1}$")
plt.yscale('log')
plt.xlim(0, 0.018)
plt.xticks([0.000, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018],
["0", "2", "4", "6", "8", "10", "12", "14", "16", "18"])
#plt.ylim(2000, 8000)
#plt.yticks([2000, 3000, 4000, 5000, 6000, 7000, 8000],
#["2", "3", "4", "5", "6", "7", "8"])
plt.errorbar(CuD, (CuN/CuT - Nulleffekt), np.sqrt(CuN/CuT - Nulleffekt),  fmt='ro', label="Messwerte Kupfer")
plt.plot(x1plot, np.exp(f(x1plot, *paramsCu)), 'b', label="Regression")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("CuGamma.pdf")

# Zink
x2plot = np.linspace(0, 21*10**(-3))
plt.figure(2)
plt.xlabel(r"$D/\mathrm{mm}$")
plt.ylabel(r"$N / \mathrm{s}^{-1}$")
plt.xlim(0, 21*10**(-3))
plt.xticks([0.000, 0.005, 0.01, 0.015, 0.02], ["0", "5", "10", "15", "20"])
#plt.ylim(2000, 7000)
#plt.yticks([2000, 3000, 4000, 5000, 6000, 7000],
#["2", "3", "4", "5", "6", "7"])
plt.yscale('log')
plt.errorbar(ZnD, (ZnN/ZnT - Nulleffekt), np.sqrt(ZnN/ZnT - Nulleffekt), fmt='ro', label="Messwerte Zink")
plt.plot(x2plot, np.exp(f(x2plot, *paramsZn)), 'b', label="Regression")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("ZnGamma.pdf")

# Auswertung b)
NL = 2.69*10**(19) # in cm^(-3)
sigma = 2.565592506*10**(-25) # in cm²

zZn = 30
rhoZN = 7.14 # g/cm^3
MZn = 65.38 # g/mol

zCu = 29
rhoCu = 8.96 # g/cm³
MCu = 63.546 # g/mol

muZn = zZn*const.N_A*rhoZN/MZn * sigma *100
muCu = zCu*const.N_A*rhoCu/MCu * sigma *10**(2)

print("")
print("Aufgabenteil b)")
print("Errechneter Wert mu Zink: ", muZn)
print("Errechneter Wert mu Kupfer: ", muCu)
print("Relative Abweichung Zink: ", (1- ZnMu/muZn)*100)

# Aufgabenteil c)
# Einlesen der Daten
AlD, AlDerr, AlT, AlN = np.genfromtxt('ALBeta.txt', unpack=True)
# Nulleffekt
NulleffektAl = 0.256
AlD = AlD *10**(-6)
ln_Al = np.log(AlN/AlT - NulleffektAl)
# Lineare Regression
def g(x, m, n):
    return m*x + n
paramsAl1, covarianceAl1 = curve_fit(g, AlD[0:4], ln_Al[0:4])
errorsAl1 = np.sqrt(np.diag(covarianceAl1))
AlN1 = ufloat(paramsAl1[1], errorsAl1[1])
#AlN1 = unp.exp(AlN1)
AlM1 = ufloat(paramsAl1[0], errorsAl1[0])

#def h(x, b):
#    return b
#paramsAl2, covarianceAl2 = curve_fit(h, AlD[8], ln_Al[8])
#errorsAl2 = np.sqrt(np.diag(covarianceAl2))
#AlN2 = ufloat(paramsAl2[0], errorsAl2[0])
#AlN2 = unp.exp(AlN2)
#AlM2 = ufloat(paramsAl2[0], errorsAl2[0])
print("")
print("Aufgabenteil c)")
print("Y-Achsenabschnitt: ",AlN1)
print("Steigung: ", AlM1)
#print(AlN2)

# Bestimmung von Rmax

r = (NulleffektAl - AlN1)/(AlM1)*10**(6)
print("Rmax in mikro m: ",r)
rho = 2.71
rmax = r*10**(-4) * rho
E = 1.92*unp.sqrt((rmax)**2 + 0.22*rmax)
print("Emax in MeV: ", E)

# Plot
plt.figure(3)
al2 = np.linspace(90*10**(-6), 277*10**(-6))
al = np.linspace(100*10**(-6), 410*10**(-6))
plt.yscale('log')
plt.xlabel(r"$D / \mathrm{\mu m}$")
plt.ylabel(r"$N / \mathrm{s}^{-1}$")
plt.xlim(90*10**(-6), 410*10**(-6))
plt.xticks([1*10**(-4), 1.5*10**(-4), 2*10**(-4), 2.5*10**(-4), 3*10**(-4), 3.5*10**(-4), 4*10**(-4)],
["100", "150", "200", "250", "300", "350", "400"])
plt.plot(AlD, (AlN/AlT), 'ro', label="Messwerte")
plt.axhline(y=NulleffektAl,xmin=0, xmax=1, hold=None, label="Untergrund")
plt.plot(al2, np.exp(g(al2, *paramsAl1)), 'g', label="Regression 1")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("AlBeta.pdf")
