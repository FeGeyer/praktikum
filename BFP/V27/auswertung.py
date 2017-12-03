import numpy as np
import matplotlib.pyplot as plt
# import uncertainties.unumpy as unp
from uncertainties import ufloat
# from scipy.optimize import curve_fit
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

# Einlesen der aus den Bildern gewonnenen Daten
ds_Rot_0, dels_Rot_0 = np.genfromtxt('Rot0.txt', unpack=True)
ds_Blau_0, dels_Blau_0 = np.genfromtxt('blau0.txt', unpack=True)
ds_Blau_90, dels_Blau_90 = np.genfromtxt('blau90.txt', unpack=True)

I_Rot_0 = 10
I_Blau_0 = 5.5
I_Blau_90 = 14

Disp_Geb_Bl = 26.95 * 10**(-12)
Disp_Geb_Ro = 48.91 * 10**(-12)


# Umrechnen in Wellenlängenänderung
def del_lambda(del_s, delta_s, Disp):
    return 0.5*(del_s/delta_s)*Disp


def mg(dl, B, lam):
    return const.h * const.c / (lam**2*B*const.value("Bohr magneton")) * dl


def B(I):
    return (x_3*I**3+x_2*I**2+x_1*I+x_0)*10**(-3)


lam_Rot = del_lambda(dels_Rot_0, ds_Rot_0, Disp_Geb_Ro)
lam_Rot_mean = ufloat(np.mean(lam_Rot), stats.sem(lam_Rot))
lam_Blau0 = del_lambda(dels_Blau_0, ds_Blau_0, Disp_Geb_Bl)
lam_Blau0_mean = ufloat(np.mean(lam_Blau0), stats.sem(lam_Blau0))
lam_Blau90 = del_lambda(dels_Blau_90, ds_Blau_90, Disp_Geb_Bl)
lam_Blau90_mean = ufloat(np.mean(lam_Blau90), stats.sem(lam_Blau90))

print(lam_Rot_mean, lam_Blau0_mean, lam_Blau90_mean)

mg_Rot = mg(lam_Rot_mean, B(I_Rot_0), 644*10**(-9))
mg_blau_0 = mg(lam_Blau0_mean, B(I_Blau_0), 480*10**(-9))
mg_blau_90 = mg(lam_Blau90_mean, B(I_Blau_90), 480*10**(-9))
print(mg_Rot, mg_blau_0, mg_blau_90)
print(B(I_Rot_0), B(I_Blau_0), B(I_Blau_90))

np.savetxt('TabelleRot.txt', np.column_stack([ds_Rot_0, dels_Rot_0, lam_Rot*10**12]), fmt="%.2f")
np.savetxt('TabelleBlau0.txt', np.column_stack([ds_Blau_0, dels_Blau_0, lam_Blau0*10**12]), fmt="%.2f")
np.savetxt('TabelleBlau90.txt', np.column_stack([ds_Blau_90, dels_Blau_90, lam_Blau90*10**12]), fmt="%.2f")
