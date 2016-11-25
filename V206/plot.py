import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 10

#Deklarieren verschiedern Veriablen. Generell gilt: 1 := wärmer werdendes Gefäß,
# T2 := kälter werdendes Gefäß
tmin, T1, Pb, T2, Pa, N = np.genfromtxt('data.txt', unpack=True)
ts = tmin*60
T1 = T1+273.15
T2 = T2+273.15
zeiten = np.array([240, 480, 720, 960])
cmapperat = 660
cwasser = 4.19e3
cmwasser = cwasser * 3

#Definition der Funktionen für Ausgleichsrechnung und Differntialquotienten
def f(x, a, b, c):
    return a * (x**2) + b * x + c

def fstrich(x, a, b):
    return 2 * a * x + b

#Ausgleichsrechnung für Temperaturverläufe und speichern der Koeefizienten A, B, und C (5b) )
params1, covariance1 = curve_fit(f, ts, T1)
errors1 = np.sqrt(np.diag(covariance1))
a1 = ufloat(params1[0], errors1[0])
b1 = ufloat(params1[1], errors1[1])
c1 = ufloat(params1[2], errors1[2])

params2, covariance2 = curve_fit(f, ts, T2)
errors2 = np.sqrt(np.diag(covariance2))
a2 = ufloat(params2[0], errors2[0])
b2 = ufloat(params2[1], errors2[1])
c2 = ufloat(params2[2], errors2[2])

#Berechnung der Quotienten (5b) )
Quotienten1 = fstrich(zeiten, a1, b1)
Quotienten2 = fstrich(zeiten, a2, b2)

#Berechnung der Güteziffern (5c) )
Temperatur1 = f(zeiten, a1, b1, c1)
Temperatur2 = f(zeiten, a2, b2, c2)
nuemp = (cmapperat + cmwasser) * Quotienten1 * (1/np.mean(N))
nuid = Temperatur1 / (Temperatur1 - Temperatur2)

#Ausgleichsrechnung zur Bestimmung der Verdampfungswärme L (5d) )
def DDK(x, L, K):
    return L * x + K
Dampfdruck = np.log(Pb)
Lparams, Lcovariance = curve_fit(DDK, Dampfdruck, T1)
Lerrors = np.sqrt(np.diag(Lcovariance))
g = ufloat(Lparams[0], Lerrors[0])
h = ufloat(Lparams[1], Lerrors[1])
print(g)
print(h)

#Ausgaben
print("T1")
print('a =', params1[0], '±', errors1[0])
print('b =', params1[1], '±', errors1[1])
print('c =', params1[2], '±', errors1[2])
print("T2")
print('a =', params2[0], '±', errors2[0])
print('b =', params2[1], '±', errors2[1])
print('c =', params2[2], '±', errors2[2])

#Plot
plt.title("Temperaturverlauf in den Gefäßen")

x_plot = np.linspace(0, 18*60)

plt.ylabel('Temperatur (K)')
plt.xlabel("Zeit (s)")
plt.plot(ts, T1, 'r+', label="T1")
plt.plot(ts, T2, 'y+', label="T2")
plt.plot(x_plot, f(x_plot, *params1), 'b-', label='Regression', linewidth=3)
plt.plot(x_plot, f(x_plot, *params2), 'b-', label='Regression', linewidth=3)
plt.legend(loc="best")
plt.savefig('Tverl.pdf')

plt.tight_layout()
