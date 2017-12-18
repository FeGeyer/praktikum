# Header
import numpy as np
import matplotlib.pyplot as plt
# import uncertainties.unumpy as unp
from uncertainties import ufloat
# from scipy.optimize import curve_fit
from scipy.stats import stats
import scipy.constants as const
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

# Daten einlesen
R_p1, R_p2, R_z1, R_z2, U, I, dt = np.genfromtxt('Messwerte.txt', unpack=True)


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


def cp(dT, U, I, dt):
    return (U * I * dt * M)/(dT * m)


dT = T_p2 - T_p1
cp = cp(dT, U, I, dt)

print(cp)
