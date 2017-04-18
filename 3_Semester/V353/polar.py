import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp

f, amp, a = np.genfromtxt('messwerte_freq.txt', unpack=True)
t, Uc = np.genfromtxt('messwerte.txt', unpack=True)
a= a*1e-6
b = 1. / f
A = amp / 20.77
U = Uc 
phi = (a / b) * 360
Phasendif = phi  * 2 * np.pi / 360

#-----------------------------------------------
def fit_Spannung(t, U):
    return -1./U * t

params, cov = curve_fit(fit_Spannung, t, np.log(U))
errors = np.diag(cov)

x = np.linspace(0, 23)
plt.plot(t, np.log(U), 'rx', label='Messdaten')
plt.plot(x, fit_Spannung(x, params), 'b-', label='Fit')
plt.ylabel(r'$\mathrm{ln\left(\frac{U_C}{U_0} \right)}$')
plt.xlabel(r'$\mathrm{t}$ / ms')
# plt.yscale('log')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('buildZeit_Amp.pdf')

print('Parameter der Zeit_Amp')
print('		n=', params[0], '\pm', errors[0])

plt.clf()
#-----------------------------------------------
def fit_Amplitude(f, RC):
    return 1. / np.sqrt(1 + f**2 * RC**2)

params, cov = curve_fit(fit_Amplitude, f, A)
errors = np.sqrt(np.diag(cov))

x = np.linspace(10, 10000)
plt.plot(f, A, 'rx', label='Messdaten')
plt.plot(x, fit_Amplitude(x, params[0]), 'b-', label='Fit')
plt.legend(loc='best')
#plt.yscale('log')
plt.xlabel(r'$\nu$ / $\mathrm{Hz}$')
plt.ylabel(r'$U_c / U_0$')
plt.savefig('buildPhase.pdf')

print('Parameter der Phase')
print('		RC=', params[0], '\pm', np.sqrt(errors[0]))

plt.clf()
#-----------------------------------------------
def fit_Phase(f, RC):
    return np.arctan(-f * RC) *360 / 2 / np.pi

params, cov = curve_fit(fit_Phase, f, phi)
errors = np.diag(cov)

x = np.logspace(0, 4)
plt.plot(f, phi, 'rx', label='Messdaten')
plt.plot(x, fit_Phase(x, params[0]), 'b-', label='Fit')
plt.xlabel(r'$\nu$ / Hz')
plt.ylabel(r'$\phi / °$')
plt.xscale('log')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('buildFrequenz_Dif.pdf')

print('Parameter der Frequenz_Dif')
print('		RC=', params[0], '\pm', errors[0])

plt.clf()

#-----------------------------------------------
def H(x, c):
    return -np.sin(c)/x

plt.polar(phi, A, 'rx')
x = np.linspace(0.0000001, 1)
plt.polar(x, H(-np.tan(x), x), 'b-')
plt.savefig('buildpolar.pdf')
plt.clf()
