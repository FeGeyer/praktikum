import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 10

#Einlesen der Daten
RbSkala, IbmA, UbV = np.genfromtxt('datab.txt', unpack=True)
RcSkala, IcmA, UcV = np.genfromtxt('datac.txt', unpack=True)
RdSkala, IRechtckmA, URechteckV, ISinusmA, USinusV = np.genfromtxt('datad.txt', unpack=True)
U0MonoGemessen = 1.55
RMessgerät = 10e6
UGegen = 3.5
#IbmA = IbmA * 10**-3

np.savetxt('B.txt', np.column_stack([((UbV/IbmA)*10**3), UbV, IbmA]),fmt="%.2f")
np.savetxt('C.txt', np.column_stack([((UcV/IcmA)*10**3), UcV, IcmA]),fmt="%.2f")
np.savetxt('Re.txt', np.column_stack([((URechteckV/IRechtckmA)*10**3), URechteckV, IRechtckmA]),fmt="%.2f")
np.savetxt('Si.txt', np.column_stack([((USinusV/ISinusmA)*10**3), USinusV, ISinusmA]),fmt="%.2f")

#Lineare Regression
def f(x, m, n):
    return m*x+n
paramsMono, covarianceMono = curve_fit(f, IbmA, UbV)
errorsMono = np.sqrt(np.diag(covarianceMono))
U0Mono = ufloat(paramsMono[1], errorsMono[1])
RiMono = ufloat(paramsMono[0], errorsMono[0])

paramsMonoGegen, covarianceMonoGegen = curve_fit(f, IcmA, UcV)
errorsMonoGegen = np.sqrt(np.diag(covarianceMonoGegen))
U0MonoGegen = ufloat(paramsMonoGegen[1], errorsMonoGegen[1])
RiMonoGegen = ufloat(paramsMonoGegen[0], errorsMonoGegen[0])

paramsRechteck, covarianceRechteck = curve_fit(f, IRechtckmA, URechteckV)
errorsRechteck = np.sqrt(np.diag(covarianceRechteck))
U0Rechteck = ufloat(paramsRechteck[1], errorsRechteck[1])
RiRechteck = ufloat(paramsRechteck[0], errorsRechteck[0])

paramsSinus, covarianceSinus = curve_fit(f, ISinusmA, USinusV)
errorsSinus = np.sqrt(np.diag(covarianceSinus))
U0Sinus = ufloat(paramsSinus[1], errorsSinus[1])
RiSinus = ufloat(paramsSinus[0], errorsSinus[0])

print("U0Mono: ", U0Mono)
print("RiMono: ", RiMono * 10**3)
print()
print("U0MonoGegen: ", U0MonoGegen)
print("RiMonoGegen: ", RiMonoGegen * 10**3)
print()
print("U0Rechteck: ", U0Rechteck)
print("RiRechteck: ", RiRechteck * 10**3)
print()
print("U0Sinus: ", U0Sinus)
print("RiSinus: ", RiSinus * 10**3)
print()
#Plotten der Klemmspannung als Funkion der Stromstärke
x1plot = np.linspace(19, 65)
plt.figure(1)
plt.title(r"Klemmspannung der Monozelle bei Gegenspannung als Funktion des Belastungsstromes")
plt.xlabel("$I/$mA")
plt.ylabel("$U_k/$V")
plt.ylim(1.14,1.41)
plt.xlim(19,65)
plt.plot(IbmA, UbV, 'b+', label='Messwerte')
plt.plot(x1plot, f(x1plot, *paramsMono) , 'g-', label='Regression')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Mono.pdf")

x2plot = np.linspace(89, 241)
plt.figure(2)
plt.title(r"Klemmspannung der Monozelle als Funktion des Belastungsstromes")
plt.xlabel("$I/$mA")
plt.ylabel("$U_k/$V")
plt.ylim(1.74,2.61)
plt.xlim(89,241)
plt.plot(IcmA, UcV, 'b+', label='Messwerte')
plt.plot(x2plot, f(x2plot, *paramsMonoGegen) , 'g-', label='Regression')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("MonoGegen.pdf")

x3plot = np.linspace(1.8, 6.9)
plt.figure(3)
plt.title(r"Klemmspannung des 1-V-Rechteckgenerators als Funktion des Belastungsstromes")
plt.xlabel("$I/$mA")
plt.ylabel("$U_k/$V")
plt.ylim(0.19,0.56)
plt.xlim(1.8, 6.9)
plt.plot(IRechtckmA, URechteckV, 'b+', label='Messwerte')
plt.plot(x3plot, f(x3plot, *paramsRechteck) , 'g-', label='Regression')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Rechteck.pdf")

x4plot = np.linspace(0.19,1.21)
plt.figure(4)
plt.title(r"Klemmspannung des 1-V-Sinusgenerators als Funktion des Belastungsstromes")
plt.xlabel("$I/$mA")
plt.ylabel("$U_k/$V")
plt.ylim(0.15, 1.25)
plt.xlim(0.19, 1.21)
plt.plot(ISinusmA, USinusV, 'b+', label='Messwerte')
plt.plot(x4plot, f(x4plot, *paramsSinus) , 'g-', label='Regression')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Sinus.pdf')

#Systematischer Fehler durch endlichen Widerstand des Voltmeters
Ri = unp.nominal_values(RiMono) * -1 * 10**3
print("Systematischer Fehler durch endl. Voltmeterwiderstand: ", (RiMono* -1 * 10**3)/RMessgerät)
print("deltaU: ", U0MonoGemessen*(RiMono* -1 * 10**3)/RMessgerät)

# Plotten der umgesetzten Leistung am Widerstand an der Monozelle
N = UbV * IbmA
Ra = UbV/(IbmA) * 10**3
NTheorie = ((Ra * U0MonoGemessen ** 2)/(Ra + Ri)**2)*10**3
np.savetxt('N.txt', np.column_stack([Ra, N, NTheorie]),fmt="%.2f")
print(Ra)
print(N)
NTheorieA = ((Ra * U0MonoGemessen ** 2)/(Ra + RiMono * -1 * 10**3)**2)*10**3
print(NTheorieA)

#print(NTheorie)
plt.figure(5)
plt.title(r"Umgesetze Leistung am Belastungswiderstand")
plt.ylabel("$N$/mW")
plt.xlabel("$R_a$/$\Omega$")
plt.xlim(17.5, 70.5)
plt.ylim(27.5, 74.1)
plt.plot(Ra, N, 'b+', label='Werte')
plt.plot(Ra, unp.nominal_values(NTheorie) , 'g-', label='Theorie')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('Leistung.pdf')
