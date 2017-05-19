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
iEinzel = iEinzel*10**(3)
iEinzel = iEinzel - 2
deltaEinzel = dEinzel - 25
print(iEinzel)
# Konstanten definieren
L = (145-20)*10
wellenlaenge = 635*10**(-6)
#print(wellenlaenge)
#print(L)
#Fit
#def f (deltax, A0, b):
#    return A0**2*b**2*np.sinc((np.pi*b*np.sin(deltax/L))/wellenlaenge)**2

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

xEinzel = np.linspace(5, 45)
plt.figure(1)
plt.xlabel(r"$\Delta x/L$")
plt.ylabel(r"$I$")
plt.plot(dEinzel, iEinzel, 'r+', label="Messwerte")
plt.plot(xEinzel, g(xEinzel, *paramsEinzel), 'b', label="Regression")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("einzel.pdf")
