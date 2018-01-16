# Header
import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
import scipy.constants as const
import scipy.integrate as integrate
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

# Daten einlesen
dt, N = np.genfromtxt('verz.txt', unpack=True)

plt.figure(1)
plt.ylabel(r"$\theta_D / T$")
plt.xlabel(r"$c_V$ / J$\,$Mol$^{-1}\,$K$^{-1}$")
plt.plot(dt, N, 'rx', label="Tabellierte Werte")
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Plateau.pdf")
plt.clf()
