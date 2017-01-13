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
finkHz, dtinmus, UcindV, UerinV  = np.genfromtxt('messwerte.txt', unpack=True)
UtiindV, tiinmus = np.genfromtxt('messwerte_obere_einh.txt', unpack=True)
dt = dtinmus * 10**(-6)
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
plt.xlabel(r"$\nu / 10^3 \, \mathrm{Hz}$")
plt.ylabel(r"$U / \mathrm{V}$")
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
plt.xlabel(r"$\nu / 10^3 \, \mathrm{Hz}$")
plt.ylabel(r"$U / \mathrm{V}$")
plt.plot(finkHz, U, 'b+', label='Messwerte')
plt.legend(loc="best")
plt.tight_layout
plt.savefig('lin.pdf')

#Teil d)
T=1/(finkHz*10**3)
phi = 360*dt/T
print(phi)

plt.figure(4)
plt.xlabel(r"$\nu / 10^3 \, \mathrm{Hz}$")
plt.ylabel(r"$\varphi / \mathrm{rad}$")
#plt.title("Spannung als Funktion des Abstandes zwischen LED und PD")
#plt.yscale("log")
plt.ylim(0,200)
plt.xlim(20,41)
plt.xscale("log")
plt.yticks([0, 45, 90, 135, 180],
           [r"$0$", r"$\frac{1}{4}\pi$", r"$\frac{1}{2}\pi$", r"$\frac{3}{4}\pi$", r"$\pi$"])
plt.plot(finkHz, phi, 'b+', label='Messwerte')
plt.legend(loc="best")
plt.tight_layout
plt.savefig('Phasehalblog.pdf')

plt.figure(5)
plt.xlabel(r"$\nu / 10^3 \, \mathrm{Hz}$")
plt.ylabel(r"$\varphi / \mathrm{rad}$")
#plt.title("Spannung als Funktion des Abstandes zwischen LED und PD")
#plt.yscale("log")
#plt.xscale("log")
plt.ylim(0,200)
plt.xlim(20,41)
plt.yticks([0, 45, 90, 135, 180],
           [r"$0$", r"$\frac{1}{4}\pi$", r"$\frac{1}{2}\pi$", r"$\frac{3}{4}\pi$", r"$\pi$"])
plt.plot(finkHz, phi, 'b+', label='Messwerte')
plt.legend(loc="best")
plt.tight_layout
plt.savefig('Phaselin.pdf')
