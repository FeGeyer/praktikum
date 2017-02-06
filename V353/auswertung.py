import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
#from uncertainties import uarray
from scipy.optimize import curve_fit
from scipy.stats import stats
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 20

#Einlesen der Daten
t_i_in_mys, U_c_von_ti_in_V = np.genfromtxt('messwerte.txt', unpack=True)
nu_in_Hz, U_c_in_V, dT_in_ms = np.genfromtxt('messwerte_freq.txt', unpack=True)

t_i_in_s = t_i_in_mys * 10 ** (-6)

#Ausgleichsrechnung a)
def f(x, a, b):
    return a * x + b

params1, covariance1 = curve_fit(f, t_i_in_s, np.log(U_c_von_ti_in_V))
errors1 = np.sqrt(np.diag(covariance1))
a = ufloat(params1[0], errors1[0])
b = ufloat(params1[1], errors1[1])

RC_a = -1 * 1/a
print(a)
print(b)
print(RC_a)
print(unp.exp(b))

U_0 = np.exp(unp.nominal_values(b))
RC = unp.nominal_values(RC_a)
print(a)
print(U_0)


plt.figure(0)
#plt.title("Frequenzabhängigkeit der Amplitude")
t_plot = np.linspace(-0.00005, 3.55*10**(-3))
plt.xlim(-0.00005, 0.00355)
plt.xlabel(r"$t_i / \mathrm{ms}$")
plt.ylabel(r"$U_{\mathrm{C}}(t_i) / \mathrm{V}$")
#plt.grid(True,which="both",ls="-")
#plt.yscale("log")
plt.xticks([0, 0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035], [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5])
plt.plot(t_i_in_mys*10**(-6), U_c_von_ti_in_V, 'r+', label="Messwerte")
plt.plot(t_plot, U_0*np.exp(-1*t_plot*1/RC), 'b-', label='Regression')
plt.legend(loc="best")
plt.tight_layout
plt.savefig('Entladung.pdf')

#Plotten b) und Ausgelichsrechnung

def g(n, a):
    return 1/np.sqrt(1 + (n * (2 * np.pi) )**2 * a**2)

params2, covariance2 = curve_fit(g, nu_in_Hz, U_c_in_V/20.77)
errors2 = np.sqrt(np.diag(covariance2))
a_b = ufloat(params2[0], errors2[0])

print(a_b)

x_plot = np.linspace(0,10000, 10000)

plt.figure(1)
#plt.title("Frequenzabhängigkeit der Amplitude")

plt.xlabel(r"$\nu / \mathrm{Hz}$")
plt.ylabel(r"$\frac{A(\nu)}{U_0}$")

plt.xscale("log")
#plt.grid(True,which="both",ls="-")
plt.plot(nu_in_Hz, U_c_in_V/20.77, 'r+', label="Messwerte")
plt.plot(x_plot, g(x_plot, *params2), 'b-', label='Regression')
plt.legend(loc="best")
plt.tight_layout
plt.savefig('Amplitude.pdf')

np.savetxt('Amplituden.txt', np.column_stack([nu_in_Hz, U_c_in_V]),fmt="%.3f")
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

plt.figure(2)
#plt.title("Frequenzabhängigkeit der Phase")

plt.xlabel(r"$\nu / \mathrm{Hz}$")
plt.ylabel(r"$\varphi / \mathrm{deg}$")

plt.xscale("log")
#plt.grid(True,which="both",ls="-")
plt.plot(nu_in_Hz, phi, 'r+', label="Messwerte")
plt.plot(x_plot, h(x_plot, *params3), 'b-', label='Regression')
plt.legend(loc="best")
plt.tight_layout
plt.savefig('Phase.pdf')

np.savetxt('Phasen.txt', np.column_stack([nu_in_Hz, dT_in_ms, phi]),fmt="%.3f")
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
