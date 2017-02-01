import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
#from uncertainties import uarray
from scipy.optimize import curve_fit
from scipy.stats import stats
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 15

#Einlesen der Daten
t_i_in_mys, U_c_von_ti_in_V = np.genfromtxt('messwerte.txt', unpack=True)
nu_in_Hz, U_c_in_V, dT_in_ms = np.genfromtxt('messwerte_freq.txt', unpack=True)

t_i_in_s = t_i_in_mys * 10 ** (-6)

#Ausgleichsrechnung a)
def f(x, a, b):
    return a * x + b

params1, covariance1 = curve_fit(f, t_i_in_s, np.log(U_c_von_ti_in_V/20.77))
errors1 = np.sqrt(np.diag(covariance1))
a = ufloat(params1[0], errors1[0])
b = ufloat(params1[1], errors1[1])

RC_a = -1 * 1/a

print(RC_a)

#Plotten b) und Ausgelichsrechnung

def g(n, a):
    return 1/np.sqrt(1 + (n * (2 * np.pi) )**2 * a**2)

params2, covariance2 = curve_fit(g, nu_in_Hz, U_c_in_V/20.77)
errors2 = np.sqrt(np.diag(covariance2))
a_b = ufloat(params2[0], errors2[0])

print(a_b)

x_plot = np.linspace(0,10000, 10000)

plt.figure(0)
plt.title("Frequenzabhängigkeit der Amplitude")

plt.xlabel(r"$\nu / \mathrm{Hz}$")
plt.ylabel(r"$\frac{A(\nu)}{U_0}$")

plt.xscale("log")
#plt.grid(True,which="both",ls="-")
plt.plot(nu_in_Hz, U_c_in_V/20.77, 'r+', label="Messwerte")
plt.plot(x_plot, g(x_plot, *params2), 'b-', label='Regression')
plt.legend(loc="best")
plt.tight_layout
plt.savefig('Amplitude.pdf')

#Plotten c) und Ausgleichsrechnung


f = 1/nu_in_Hz
phi = dT_in_ms * 10 ** (-3) * 360 / f

def h(n, a):
    return np.arctan(-1 * n * (2*np.pi) * a) * 360 / np.pi / 2

params3, covariance3 = curve_fit(h, nu_in_Hz, phi)
errors3 = np.sqrt(np.diag(covariance3))
a_c = ufloat(params3[0], errors3[0])

print(a_c)
print(phi)

plt.figure(1)
plt.title("Frequenzabhängigkeit der Phase")

plt.xlabel(r"$\nu / \mathrm{Hz}$")
plt.ylabel(r"$\varphi$")

plt.xscale("log")
#plt.grid(True,which="both",ls="-")
plt.plot(nu_in_Hz, phi, 'r+', label="Messwerte")
plt.plot(x_plot, h(x_plot, *params3), 'b-', label='Regression')
plt.legend(loc="best")
plt.tight_layout
plt.savefig('Phase.pdf')

#Polarplot
def H(x, c):
    return -np.sin(c)/x

phi = dT_in_ms * 10 ** (-3) * 2 * np.pi * nu_in_Hz
A = U_c_in_V/20.77
plt.figure(2)
plt.polar(phi, A, 'rx', label="Messwerte")
x = np.linspace(0.0000001,  2, 1000)
plt.polar(x, H(-np.tan(x), x), 'b-', label="Theoriekurve")
plt.legend(loc="lower left")
plt.tight_layout
plt.savefig('Polar.pdf')
