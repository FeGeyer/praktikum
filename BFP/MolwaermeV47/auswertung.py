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

# Molwärmen ausrechnen, Molmasse berechnen.
u = const.value("atomic mass constant")
M = ufloat(63.546, 0.003) * u * const.Avogadro
print(M)
m = 0.342


# cp berechnen
def cp(dT, U, I, M, dt, m):
    return (U * I * dt * M)/(dT * m)


dT = T_p2 - T_p1
cp = cp(dT, U, I, M, dt, m)

np.savetxt('TexTabellen/cp.txt', np.column_stack([
        unp.nominal_values(T_p1)+273.15,
        unp.std_devs(T_p1),
        unp.nominal_values(T_p2)+273.15,
        unp.std_devs(T_p2),
        unp.nominal_values(dT),
        unp.std_devs(dT),
        U,
        I * 10**3,
        unp.nominal_values(dt),
        unp.std_devs(dt),
        unp.nominal_values(cp),
        unp.std_devs(cp)
        ]), delimiter=' & ', newline=r' \\'+'\n', fmt='%.2f')


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
plt.plot(T_Koeff, Koeff * 10**6, 'r^', label="Tabellierte Werte")
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
V0 = 7.10 * 10 ** (-6)

cv = cv(cp, alpha_T, kappa, V0, Tbar)

np.savetxt('TexTabellen/cv.txt', np.column_stack([
        unp.nominal_values(Tbar),
        unp.std_devs(Tbar),
        unp.nominal_values(alpha_T) * 10**(6),
        unp.std_devs(alpha_T) * 10**(6),
        unp.nominal_values(cp),
        unp.std_devs(cp),
        unp.nominal_values(cv),
        unp.std_devs(cv)
        ]), delimiter=' & ', newline=r' \\'+'\n', fmt='%.2f')

plt.figure(2)
plt.xlim(unp.nominal_values(Tbar).min()-2, unp.nominal_values(Tbar).max()+2)
plt.ylim(unp.nominal_values(cv).min()-1, unp.nominal_values(cp).max()+1)
plt.ylabel(r"$c$ / J$\,$Mol$^{-1}\,$K$^{-1}$")
plt.xlabel(r"$T$ / K")

# plt.errorbar(x=unp.nominal_values(Tbar), y=unp.nominal_values(cp),
#              xerr=unp.std_devs(Tbar),
#              yerr=unp.std_devs(cv), fmt='b^', label=R"$c_\mathrm{p}$")
plt.plot(unp.nominal_values(Tbar), unp.nominal_values(cp), 'bv',
         label=r"$c_\mathrm{p}$")
# plt.plot(unp.nominal_values(Tbar), unp.nominal_values(cp), 'b--', )

# plt.errorbar(x=unp.nominal_values(Tbar), y=unp.nominal_values(cv),
#              xerr=unp.std_devs(Tbar),
#              yerr=unp.std_devs(cv), fmt='ro', label=r"$c_\mathrm{V}$")
plt.plot(unp.nominal_values(Tbar), unp.nominal_values(cv), 'r^',
         label=r"$c_\mathrm{V}$")
# plt.plot(unp.nominal_values(Tbar), unp.nominal_values(cv), 'r--', )

plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Kapaz.pdf")
plt.clf()

# Theta_D / T Werte für die c_V Werte bestimmen. Hierbei werden der für unsere
# Daten relevante Bereich aus der Tabelle abgelesen und entsprechend gefittet.

# Es sollen nur Werte bis 170K betrachtet werden
cv_bis_170 = cv[Tbar <= 170]
Tbar_bis_170 = Tbar[Tbar <= 170]
print(cv_bis_170)

# Einlesen, Fitten, zur Überprüfung plotten
theta_T_Tabelle, c_V_Tabelle = np.genfromtxt('debeye.txt', unpack=True)
x2 = np.linspace(13, 20)

debeye, cov2 = np.polyfit(c_V_Tabelle, theta_T_Tabelle, 1, cov=True, )
errors2 = np.sqrt(np.diag(cov2))
d_1 = ufloat(debeye[0], errors[0])
d_0 = ufloat(debeye[1], errors[1])
Fit2 = np.polyval(debeye, x2)

print(debeye)
print(d_0)
print(d_1)

plt.figure(3)
plt.ylim(2.0, 3.9)
plt.xlim(13, 20)
plt.ylabel(r"$\theta_D / T$")
plt.xlabel(r"$c_V$ / J$\,$Mol$^{-1}\,$K$^{-1}$")
plt.plot(c_V_Tabelle, theta_T_Tabelle, 'r^', label="Tabellierte Werte")
plt.plot(x2, Fit2, 'r--', label="Regression")
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Debeye.pdf")
plt.clf()

# Fitfunktion auswerten
theta_d_durch_T = np.polyval(debeye, cv_bis_170)
theta_d = theta_d_durch_T * Tbar_bis_170

# Mittelwerte von theta_D:

# Zuerst: Good ol' arithmetisches Mittel
theta_d_arith = np.mean(theta_d)
print(theta_d_arith)

# Nun: Ein Gewichteter Mittelwert. Allgemein: np.sum(w * x), wobei w das array
# der Gewichte ist und x die eigentlichen werte sind. Die Summe der Gewichte
# muss 1 sein. Jetzt muss man nur noch Gewichte finden. Hier werden nun Werte
# weniger stark gewichtet, bei denen Zylinder und Probe stark unterschiedlich
# temperiert waren.

# Mittelwerte der Zylindertemperaturen:
TZylbar = T_z1 + 0.5*(T_z2 - T_z1) + 273.15
T_z1_bis_170 = T_z1[Tbar <= 170] + 273.15
T_z2_bis_170 = T_z2[Tbar <= 170] + 273.15
TZylbar_bis_170 = TZylbar[Tbar <= 170]

rel = Tbar_bis_170 / TZylbar_bis_170

# Durchloopen, um Fälle zu filtern, in denen Tbar > Tzylbar war zu finden und
# umzudrehen
for i in range(len(rel)):
    if rel[i] > 1:
        rel[i] = TZylbar_bis_170[i] / Tbar_bis_170[i]

# rel soll das Gewichtsarray werden. Also: Auf 1 normieren.
summe = np.sum(rel)
rel = rel / np.sum(rel)

# gewichteten Mittelwert ausrechnen:
theta_d_gew = np.sum(rel*theta_d)
print(theta_d_gew)

# Jetzt Teil d: Zuerst brauchen wir N_L und V:
N_L = m/M * const.Avogadro
V = V0 * m/M
v_long = 4700
v_trans = 2260

# Aus V L bestimmen
L = V**(1/3)

# In die Fomel einsetzten:
w_D = ((18*np.pi**2*N_L) / (L**3*((1/v_long**3) + (2/v_trans**3))))**(1/3)

print(w_D*10**(-12), "in THz")

# Daraus theta_D
theta_D_2 = const.hbar * w_D / const.k

print(theta_D_2)

print(N_L)
print(L)

rel = rel * summe

np.savetxt('TexTabellen/theta.txt', np.column_stack([
        unp.nominal_values(Tbar_bis_170),
        unp.std_devs(Tbar_bis_170),
        unp.nominal_values(T_z1_bis_170),
        unp.std_devs(T_z1_bis_170),
        unp.nominal_values(T_z2_bis_170),
        unp.std_devs(T_z2_bis_170),
        unp.nominal_values(TZylbar_bis_170),
        unp.std_devs(TZylbar_bis_170),
        unp.nominal_values(rel) * 100,
        unp.std_devs(rel) * 100,
        unp.nominal_values(cv_bis_170),
        unp.std_devs(cv_bis_170),
        unp.nominal_values(theta_d),
        unp.std_devs(theta_d)
        ]), delimiter=' & ', newline=r' \\'+'\n', fmt='%.2f')
