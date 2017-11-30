import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
import scipy.constants as const
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

# Magnetfeldeichung
# Einlesen der Daten
I, B = np.genfromtxt('B_Feld.txt', unpack=True)

# Fitten
x1 = np.linspace(-0.25, 20.25)
B_Koeff, cov = np.polyfit(I, B, 3, cov=True, )
errors = np.sqrt(np.diag(cov))
x_3 = ufloat(B_Koeff[0], errors[0])
x_2 = ufloat(B_Koeff[1], errors[1])
x_1 = ufloat(B_Koeff[2], errors[2])
x_0 = ufloat(B_Koeff[3], errors[3])
Fit = np.polyval(B_Koeff, x1)
B_Koeff = np.array([x_3, x_2, x_1, x_0])
print(B_Koeff)

# Plotten
plt.figure(1)
plt.xlim(-0.25, 20.25)
plt.ylim(0, 1050)
plt.xlabel(r"$I / \mathrm{A}$")
plt.ylabel(r"$B / \mathrm{mT}$")
plt.plot(I, B, 'rx', label="Messwerte")
plt.plot(x1, Fit, 'r--', label="Regression")
plt.xticks([0, 5, 10, 15, 20], ["0", "5", "10", "15", "20"])
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("B_Feld.pdf")
plt.clf()
