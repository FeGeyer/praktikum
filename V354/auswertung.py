import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 10

#Teil a)
#Einlesen der Daten
finkHz, tinmus, UcindV, UerinV  = np.genfromtxt('messwerte.txt', unpack=True)
UtiindV, tiinmus = np.genfromtxt('messwerte_obere_einh.txt', unpack=True)
Uc = UcindV/10
Uti = UtiindV/10
tiins = tiinmus * 10**(-6)
#Fitten der e-Funktion
def f(t, A, mu):
    return A*np.exp(mu*t)
params, covariance = curve_fit(f, tiins, Uti)
errors = np.sqrt(np.diag(covariance))
A = ufloat(params[0], errors[0])
mu = ufloat(params[1], errors[1])

#Berechnen von R unt T und Ausgabe
L=ufloat(10.11, 0.03)*10**(-3)
mu=mu/(-2*np.pi)
R = (mu * 4 * np.pi * L)
T = 2*L/R
print("Ao: ",A)
print("mu: ",mu)
print("Reff: ", R)
print("Teff: ",T)

#Plotten der e-Funktion
x1plot = np.linspace(4*10**(-5), 3*10**(-4))
plt.figure(1)
plt.title("Spannung als Funktion des Abstandes zwischen LED und PD")
#plt.yscale("log")
#plt.xscale("log")
plt.plot(tiinmus, Uti, 'b+', label='Gefittete Messwerte')
plt.plot(x1plot, f(x1plot, *params) , 'g-', label='Regression')
plt.legend(loc="best")
plt.tight_layout
plt.savefig('Abst.pdf')

#Teil b)
#Berechnen von Rap
C = ufloat(2.098, 0.006) * 10**(-9)
Rap = unp.sqrt((4*L)/(C))
print("Rap: ",Rap)

#Teil c)
#Halblogarithmische Darstellung der Kondensatorspannung U=Uc/Uer
U = Uc/UerinV
print(Uc)
print(UerinV)
print(U)
R2 = ufloat(509.5, 0.5)
w0 = unp.sqrt(1/(L*C))
qerr = 1/(w0*R2*C)
print(qerr)
plt.figure(2)
#plt.title("Spannung als Funktion des Abstandes zwischen LED und PD")
#plt.yscale("log")
plt.xscale("log")
plt.plot(finkHz, U, 'b+', label='Messwerte')
plt.legend(loc="best")
plt.tight_layout
plt.savefig('Halblog.pdf')

plt.figure(3)
#plt.title("Spannung als Funktion des Abstandes zwischen LED und PD")
#plt.yscale("log")
#plt.xscale("log")
plt.plot(finkHz, U, 'b+', label='Messwerte')
plt.legend(loc="best")
plt.tight_layout
plt.savefig('lin.pdf')
