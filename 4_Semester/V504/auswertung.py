import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
import scipy.constants as const
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

# Einlesen
Ukenn, I1, I2, I3, I4, I5 = np.genfromtxt('kennlinien.txt', unpack=True)

# Plot
plt.figure(1)
plt.xlim(0, 260)
plt.plot(Ukenn, I1,'.', label="Kennlinie 1")
plt.plot(Ukenn, I2,'.', label="Kennlinie 2")
plt.plot(Ukenn, I3,'.', label="Kennlinie 3")
plt.plot(Ukenn, I4,'.', label="Kennlinie 4")
plt.plot(Ukenn, I5,'.', label="Kennlinie 5")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("kennlinien.pdf")

# Sättigungsstrom
Is1 = np.max(I1)
Is2 = np.max(I2)
Is3 = np.max(I3)
Is4 = np.max(I4)
Is5 = np.max(I5)

print("")
print("Aufgabenteil a): ")
print(Is1)
print(Is2)
print(Is3)
print(Is4)
print(Is5)

# Aufgabenteil b)
# Werte logarithmieren
logU = np.log(Ukenn[1:7])
logI = np.log(I5[1:7])
def f (x, m, b):
    return m * x + b
paramsRaum, covarianceRaum = curve_fit(f, logU, logI, p0=[1.5, 2])
errorsRaum = np.sqrt(np.diag(covarianceRaum))
exponent = ufloat(paramsRaum[0], errorsRaum[0])
print("")
print("Aufgabenteil b): ")
print(exponent)


plt.figure(2)
x1 = np.linspace(2, 4.2)
plt.xlim(2, 4.2)
plt.xlabel(r"$U$")
plt.ylabel(r"$I$")
#plt.xscale('log')
#plt.yscale('log')
plt.plot(logU, logI,'r.', label="Raumladungsgebiet")
plt.plot(x1, f(x1, *paramsRaum), 'b-', label="Regression")
#plt.plot(np.log(Ukenn[1:25]), np.log(I5[1:25]), 'g.', label="Alle Werte")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("raumladung.pdf")

# Aufgabenteil c)
U, I6 = np.genfromtxt('anlaufstrom.txt', unpack=True)
I6 = I6 * 10**(-9) # Umrechnen in Ampere
U = U - 1*10**6 * I6
logI6 = np.log(I6)
U = U[0:9]
logI6 = logI6[0:9]

paramsA, covarianceA = curve_fit(f, U, logI6)
errorsA= np.sqrt(np.diag(covarianceA))
exponent2 = ufloat(paramsA[0], errorsA[0])

T = - (const.e)/(const.k*exponent2)
print("")
print("Aufgabenteil c): ")
print(T)

plt.figure(3)
x2 = np.linspace(-0.1, 0.9)
plt.xlim(-0.1, 0.9)
plt.xlabel(r"$U / \mathrm{V}$")
plt.ylabel(r"$I$")
plt.plot(U, logI6, 'r.', label="Anlaufgebiet")
plt.plot(x2, f(x2, *paramsA), 'b-', label="Regression")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("anlauf.pdf")
