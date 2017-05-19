import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
import scipy.constants as const
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

#Einlesen der Daten
#Einfach Spalt
dEinzel, iEinzel = np.genfromtxt('einzel.txt', unpack=True)
# Umrechnen in nA
iEinzel = iEinzel*10**(3)
iEinzel = iEinzel - 2
# Konstanten definieren
L = (145-20)*10
wellenlaenge = 635*10**(-6)

#Fit

def g (x, A0, b, c):
    return  (A0**2) * (b**2) * (wellenlaenge/(np.pi*b*np.sin((abs(x-c))/L)))**2 * (np.sin((np.pi*b*np.sin((abs(x-c))/L))/wellenlaenge))**2
paramsEinzel, covarianceEinzel = curve_fit(g, dEinzel, iEinzel, p0=[np.max(iEinzel), 0.075, 25])
errorsEinzel = np.sqrt(np.diag(covarianceEinzel))
b = ufloat(paramsEinzel[1], errorsEinzel[1])
A0 = ufloat(paramsEinzel[0], errorsEinzel[0])
c = ufloat(paramsEinzel[2], errorsEinzel[2])
print(b)
print(A0)
print(c)

xEinzel = np.linspace(4, 46)
plt.figure(1)
plt.xlim(4, 46)
plt.xlabel(r"$x/\mathrm{mm}$")
plt.ylabel(r"$I/\mathrm{nA}$")
plt.plot(dEinzel, iEinzel, 'r+', label="Messwerte")
plt.plot(xEinzel, g(xEinzel, *paramsEinzel), 'b', label="Regression")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("einzel.pdf")

# Erster Doppelspalt klein
# Einlesen
dDoppelk, iDoppelk = np.genfromtxt('doppel_klein.txt', unpack=True)
# Umrechnen in nA
iDoppelk =iDoppelk*10**(3)
iDoppelk = iDoppelk - 2
# Fit
def f (x, A0, b, c, g):
    return A0 * (np.cos((np.pi*g*np.sin(abs(x-c)/L))/wellenlaenge))**2 * (A0**2) * (b**2) * (wellenlaenge/(np.pi*b*np.sin((abs(x-c))/L)))**2 * (np.sin((np.pi*b*np.sin((abs(x-c))/L))/wellenlaenge))**2

paramsDoppelk, covarianceDoppelk = curve_fit(f, dDoppelk, iDoppelk, p0=[np.max(iDoppelk), 0.15, 25, 0.25])
errorsDoppelk = np.sqrt(np.diag(covarianceDoppelk))
b2 = ufloat(paramsDoppelk[1], errorsDoppelk[1])
A02 = ufloat(paramsDoppelk[0], errorsDoppelk[0])
c2 = ufloat(paramsDoppelk[2], errorsDoppelk[2])
d1 = ufloat(paramsDoppelk[3], errorsDoppelk[3])
print("")
#print(np.max(iDoppelk))
print(b2)
print(A02)
print(c2)
print(d1)


# Plot
xDoppelk = np.linspace(17, 33, 1000)
plt.figure(2)
plt.xlim(17, 33)
plt.ylim(0, 5000)
plt.yticks([0, 1000, 2000, 3000, 4000, 5000],["0", "1", "2", "3", "4", "5"])
plt.xlabel(r"$x/\mathrm{mm}$")
plt.ylabel(r"$I/\mathrm{\mu A}$")
plt.plot(dDoppelk, iDoppelk, 'r+', label="Messwerte")
plt.plot(xDoppelk, f(xDoppelk, *paramsDoppelk), 'b', label="Regression")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("doppelk.pdf")
