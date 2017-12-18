# Header
import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
# from scipy.optimize import curve_fit
from scipy.stats import stats
import scipy.constants as const
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

# Daten einlesen
R_p1, R_p2, R_z1, R_z2, U, I, dt = np.genfromtxt('Messwerte.txt', unpack=True)

# Fehler auf Werte rechnen
dt = unp.uarray(dt, 5)
R_p1 = unp.uarray(R_p1, R_p1 * 0.002)
R_p2 = unp.uarray(R_p2, R_p2 * 0.002)
R_z1 = unp.uarray(R_z1, R_z1 * 0.002)
R_z2 = unp.uarray(R_z2, R_z2 * 0.002)


# Widerstände in Temperataturen umrechnen
def T(R):
    return 0.00134 * R**2 + 2.296 * R - 243.02


T_p1 = T(R_p1)
T_p2 = T(R_p2)
T_z1 = T(R_z1)
T_z2 = T(R_z2)


# Stromstärken umrechnen
I *= 10**(-3)

# Molwärmen ausrechnen, Molmasse wird von g/mol in kg/mol umgerechnent.
M = ufloat(63.546, 0.003) * 10**(-3)
m = 0.342


# cp berechnen
def cp(dT, U, I, dt):
    return (U * I * dt * M)/(dT * m)


dT = T_p2 - T_p1
cp = cp(dT, U, I, dt)


# Mittelwert der Probentemperaturdifferenzen in Kelvin
def Mittel(T_1, T_2):
    return T_1 + 0.5*(T_2-T_1) + 273.15


Tbar = Mittel(T_p1, T_p2)

# Ausdehnungskoeffizient bestimmen
T_Koeff, Koeff = np.genfromtxt('ausdehnungskoeffizient.txt', unpack=True)
Koeff *= 10**(-6)

alpha_i, cov = np.polyfit(T_Koeff, Koeff, 4, cov=True, )
errors = np.sqrt(np.diag(cov))
x_4 = ufloat(alpha_i[0], errors[0])
x_3 = ufloat(alpha_i[1], errors[1])
x_2 = ufloat(alpha_i[2], errors[2])
x_1 = ufloat(alpha_i[3], errors[3])
x_0 = ufloat(alpha_i[4], errors[4])

x1 = np.linspace(70-1, 300+1)
Fit = np.polyval(alpha_i, x1)
alpha_i = np.array([x_4, x_3, x_2, x_1, x_0])
print(alpha_i)

# Zur Überprüfung plotten
plt.figure(1)
plt.xlim(70-1, 300+1)
plt.ylim(6.5, 17)
plt.xlabel(r"$T / \mathrm{K}$")
plt.ylabel(r"$\alpha / 10^{-6} \mathrm{grd}^{-1}$")
plt.plot(T_Koeff, Koeff * 10**6, 'rx', label="Messwerte")
plt.plot(x1, Fit * 10**6, 'r--', label="Regression")
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Koeff.pdf")
plt.clf()


# cv bestimmen
def cv(cp, alpha_T, kappa, V0, Tbar):
    return -9*alpha_T**2*kappa*V0*Tbar + cp


alpha_T = np.polyval(alpha_i, Tbar)
kappa = 137.8 * 10**(9)
V0 = 7.11 * 10 ** (-6)

cv = cv(cp, alpha_T, kappa, V0, Tbar)
print(cv)
