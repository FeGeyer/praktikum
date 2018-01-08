import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
# from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
# from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
x = np.linspace(0, 3, 1000)
me = const.electron_mass
n = 3.3543  # Brechungsindex für GaAs bei einer Wellenlänge von 1771.14
B = ufloat(231, 4)
B *= 10**3  # Umrechnen auf T
A = (const.e**3) * B / (8*(np.pi**2) * const.epsilon_0*(const.c**3))    # konst zum bestimmen der eff Masse


def mittel(x):  # the real mean()-ing of life
    return ufloat(np.mean(x), np.std(x, ddof=1)/np.sqrt(len(x)))


def relf(l, m):  # in Prozent
    return (np.absolute(l-m)/l)*100


# d = dicke= Länge im Strahlengang
hd = 5.11*10**(-3)
d2 = 1.296*10**(-3)
d1 = 1.36*10**(-3)
# Werte einelsen l = Wellenläng. tb= teta mit bfeld , to= teta ohne (b/o)bm= bogenmin mit und ohne bfeld
# hochreines GaAg
hl , htb , hbbm , hto , hobm = np.genfromtxt('hochreinesGaAs.txt', unpack = True, dtype='float64')
#1. n-dotiertes GaAs
l1 , tb1, bbm1 , to1 , obm1 = np.genfromtxt('1.n-dotiertesGaAs.txt', unpack = True, dtype='float64')
#2. n-dotiertes GaAs
l2  , tb2 , bbm2 , to2 , obm2 = np.genfromtxt('2.n-dotiertesGaAs.txt', unpack = True, dtype='float64')
#jetz noch die Winkel um und zusammenrechnen. dann die ohne bFeld von den mit BFeld abziehen. Alle Messwerte in ein Plot mit 'o-'
#Winkel zusammenrechnen eine bogenmin = (1°/60) Grad
u = 1/60 #u der kleine Umrechnungsfaktor, weil faul
htb += hbbm*u
hto += hobm*u
tb1 += bbm1*u
to1 += obm1*u
tb2 += bbm2*u
to2 += obm2*u
# The True Teta ausrechnen t mit BFeld - t ohne
ht = 0.5*(htb - hto)
t1 = 0.5*(tb1 - to1)
t2 = 0.5*(tb2 - to2)
ht = abs(ht)
t1 = abs(t1)
t2 = abs(t2)
# Umrechnen der Wellenlängen
# hl*=10**-6
# l1*=10**-6
# l2*=10**-6
# Differenzzwischen der hochreinen und den zwei Proben
D1 = abs(ht/hd-t1/d1)
D2 = abs(ht/hd-t2/d2)
D = D2[0:3]
D = np.append(D, D2[5:])
h = hl[0:3]
h = np.append(h, hl[5:])
# print(D2)
# print(D)
# Fit
N1 = A*2.8*10**16   # weil umgerechnet auf Meter
N2 = A*1.2*10**16


def fitf1(x, m, b):
    return x*m + b


params1, cov1 = curve_fit(fitf1, hl**2, D1)
errors1 = np.sqrt(np.diag(cov1))
params2, cov2 = curve_fit(fitf1, hl**2, D2)
errors2 = np.sqrt(np.diag(cov2))
m1 = ufloat(params1[0], errors1[0])
b1 = ufloat(params1[1], errors1[1])
m2 = ufloat(params2[0], errors2[0])
b2 = ufloat(params2[1], errors2[1])


# ----------------------------------
def fitn(x, a, b, c):
    return a/x**2 + b/x + c


lambda_n, n, k = np.genfromtxt('n.txt', unpack=True)
lambda_n *= 10**(-3)   # Umrechnen in nm
print("Lambda: ", lambda_n)
params_n, cov_n = curve_fit(fitn, lambda_n, n)
errors_n = np.sqrt(np.diag(cov_n))
a_n = ufloat(params_n[0], errors_n[0])
b_n = ufloat(params_n[1], errors_n[1])
c_n = ufloat(params_n[2], errors_n[2])
print(a_n)
print(b_n)
print(c_n)

x_n = np.linspace(1, 2.7)
plt.plot(lambda_n, n, 'rx', label="Literaturwerte")
plt.plot(x_n, fitn(x_n, *params_n), 'b-', label="Regression")
plt.plot(hl, fitn(hl, *params_n), 'g^', label="Wellenlängen")
plt.ylabel(r"$n$")
plt.xlabel(r"$\lambda \, / \, \mathrm{\mu m}$")
plt.grid()
plt.tight_layout()
plt.legend(loc="best")
# plt.show()
plt.savefig("dispersion.pdf")
plt.clf()

x = np.linspace(1, 7, 10000)
plt.plot(hl**2, D1, 'r+', label=r"$\theta_1$")
plt.plot(x, fitf1(x, *params1), 'r-', label="Regression 1")
plt.plot(hl**2, D2, 'b+', label=r"$\theta_2$")
plt.plot(x, fitf1(x, *params2), 'b-', label="Regression 2")
plt.xlabel(r"$\lambda^2 / \mathrm{\mu m}^2$")
plt.ylabel(r"$\theta_\mathrm{diff} / \mathrm{\frac{1}{m}}$")
plt.tight_layout()
plt.grid()
plt.legend(loc="best")
plt.savefig("Test.pdf")
plt.clf()

plt.plot(hl**2, ht, 'r+', label="hochreines GaAs")
plt.plot(hl**2, t1, 'b+', label="1. n-dotiertes GaAs")
plt.plot(hl**2, t2, 'g+', label="2. n-dotiertes GaAs")
plt.xlabel(r"$\lambda^2 / \mathrm{\mu m}^2$")
plt.ylabel(r"$\theta_\mathrm{norm} / \mathrm{\frac{1}{m}}$")
plt.tight_layout()
plt.grid()
plt.legend(loc="best")
plt.savefig("b.pdf")
plt.clf()

n_l = fitn(hl, *params_n)
print(n_l)
# print("Erste Steigung: ", m1, b1)
# print("Zweite Steigung: ", m2, b2)
# print('Hildegard1:', unp.sqrt((N1/(m1*10**12*n_l))))
# print('Hildegard2:', unp.sqrt((N2/(m2*10**12*n_l))))
# print('Hildegard1/me:', unp.sqrt((N1/(m1*10**12*n_l)))/me)
# print('Hildegard2/me:', unp.sqrt((N2/(m2*10**12*n_l)))/me)
m_eff1 = unp.sqrt((N1/(m1*10**12*n_l)))
print('Mittelwert für m_eff1: ', np.mean(m_eff1))
m_eff2 = unp.sqrt((N2/(m2*10**12*n_l)))
print('Mittelwert für m_eff2: ', np.mean(m_eff2))
lol = np.array([m_eff1/me, m_eff2/me])
print("Gesamtmittelwert: ", np.mean(lol))

# Tabelle
np.savetxt('TexTabellen/HGaAstab.txt',np.column_stack([hl, htb, hto, ht,(ht/hd)]),
           delimiter=' & ',newline= r' \\'+'\n', fmt='%.2f')

np.savetxt('TexTabellen/1.nGaAstab.txt',np.column_stack([l1, tb1, to1, t1, (t1/d1)]),
           delimiter=' & ',newline= r' \\'+'\n', fmt='%.2f')

np.savetxt('TexTabellen/2.nGaAstab.txt',np.column_stack([l2, tb2, to2, t2,( t2/d2)]),
           delimiter=' & ',newline= r' \\'+'\n', fmt='%.2f')

np.savetxt('TexTabellen/Differenz2.txt',np.column_stack([D1, D2]),
           delimiter=' & ',newline= r' \\'+'\n', fmt='%.2f')

np.savetxt('TexTabellen/na.txt', np.column_stack([hl, n_l,
           unp.nominal_values(m_eff1)*10**33, unp.std_devs(m_eff1)*10**33,
           unp.nominal_values(m_eff2)*10**33, unp.std_devs(m_eff2)*10**33]),
           delimiter=' & ',newline= r' \\'+'\n', fmt='%.2f')
