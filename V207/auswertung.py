import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 10

#Einlesen der Daten
DMK1mm, DMK2mm, GK1g, GK2g = np.genfromtxt('dataDurchuGew.txt', unpack=True)
RaumTK1, RaumTK2 = np.genfromtxt('Raumtemperatur.txt', unpack=True)
Temperatur, Messung1, Messung2 = np.genfromtxt('Temperaturabh.txt', unpack=True)
TK = Temperatur + 273.15
DMK1 = DMK1mm * 10**(-3)
DMK2 = DMK2mm * 10**(-3)
GK1 = GK1g * 10**(-3)
GK2 = GK1g * 10**(-3)

#Gewicht und Durchmesser in mm und gram
Durchm1mm = ufloat(np.mean(DMK1mm), stats.sem(DMK1mm))
Durchm2mm = ufloat(np.mean(DMK2mm), stats.sem(DMK2mm))
Gew1g = ufloat(np.mean(GK1g), stats.sem(GK1g))
Gew2g = ufloat(np.mean(GK2g), stats.sem(GK2g))
RaumTK1Mean = ufloat(np.mean(RaumTK1), stats.sem(RaumTK1))
RaumTK2Mean = ufloat(np.mean(RaumTK2), stats.sem(RaumTK2))
print(RaumTK1Mean)
print(RaumTK2Mean)
print(Durchm1mm, Durchm2mm)
print(Gew1g, Gew2g)

#a) Bestimmen der Dichten der Kugeln, K1=kleinere, leichtere Kugel
Durchm1 = ufloat(np.mean(DMK1), stats.sem(DMK1))
Durchm2 = ufloat(np.mean(DMK2), stats.sem(DMK2))
Gew1 = ufloat(np.mean(GK1), stats.sem(GK1))
Gew2 = ufloat(np.mean(GK2), stats.sem(GK2))
print(Durchm1)
VK1 = (4/3) * np.pi * ((Durchm1/2) ** 3)
VK2 = (4/3) * np.pi * ((Durchm2/2) ** 3)

np.savetxt('DurchmesseruGewicht.txt', np.column_stack([DMK1mm, GK1g, DMK2mm, GK2g]), fmt="%.2f")

rho1 = Gew1/VK1
rho2 = Gew2/VK2

print("Dichte der ersten Kugel: ", rho1, " kg*m^-3")
print("Dichte der zweite Kugel: ", rho2, " kg*m^-3")
#b) Bestimmen der Aperaturkonstante für die große Kugel
FallzeitRaum1 = ufloat(np.mean(RaumTK1), stats.sem(RaumTK1))
FallzeitRaum2 = ufloat(np.mean(RaumTK2), stats.sem(RaumTK2))
rhoH2O = 997.05
Kkl = 7.64e-8
eta = Kkl * (rho1 - rhoH2O) * FallzeitRaum1
print("Eta Klein: ", eta)
Kgr = eta / ((rho2 - rhoH2O) * FallzeitRaum2)

print("Apperaturkonstante für die große Kugel: ", Kgr, " Pa*kg*m^-3")

#c) Bestimmen der Konstanten und Plot + Fit
Messwerte = (Messung1+Messung2)/2
etaT = Kgr * (rho2 - rhoH2O) * Messwerte
lnetaT = np.log(unp.nominal_values(etaT))

def f(x, m, b):
    return m * x + b
params, covariance = curve_fit(f, 1/TK, lnetaT)
errors = np.sqrt(np.diag(covariance))
m = ufloat(params[0], errors[0])
b = ufloat(params[1], errors[1])
print("Fit m: ", m)
print("Fit b: ", b)

#Reynolszahl bei der Messung der großen Kugel
v = 0.1/Messwerte
EtaBerechnetln = m * 1/TK + b
EtaBerechnet = unp.exp(EtaBerechnetln)
Re = (rho2 * v * Durchm2)/EtaBerechnet
print("Viskosität: ", EtaBerechnet)
print("Reynoldszahl: ", Re)
print(EtaBerechnet)
np.savetxt('Messwerte.txt', np.column_stack([TK, Messung1, Messung2, Messwerte, unp.nominal_values(etaT)*10**3]),fmt="%.2f")

#c) Plot
x_plot = np.linspace(0.00291418, 0.00335402)

plt.figure(0)
plt.title("ln($\eta$) als Funktion von $1/T$")

plt.ylabel('$ln(\eta)$')
plt.xlabel("$1/T$ (1/K)")
plt.plot(1/TK, lnetaT, 'r+', label="$ln(\,\eta\,(1/T\,))$")
plt.plot(x_plot, f(x_plot, *params), 'b-', label='Regression')
plt.legend(loc="best")
plt.tight_layout
plt.savefig('etaT.pdf')
