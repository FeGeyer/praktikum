import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
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

print("Kupfer: ", CuN0, CuMu)
print("Zink: ",ZnN0, ZnMu)

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
