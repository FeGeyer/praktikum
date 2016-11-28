import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 10

#Deklarieren verschiedern Veriablen. Generell gilt: 1 := wärmer werdendes Gefäß,
# T2 := kälter werdendes Gefäß
tmin, T1, Pb, T2, Pa, N = np.genfromtxt('data.txt', unpack=True)
ts = tmin*60
T1 = T1+273.15
T2 = T2+273.15
zeiten = np.array([240, 480, 720, 960])
ADruck = np.array([3.2, 3.2, 3.2, 3.2])*100000
BDruck = np.array([8, 9.5, 11.25, 12.75])*100000
cmapperat = 660
cwasser = 4.19e3
cmwasser = cwasser * 3
rhogasnorm = 5.51
kappa = 1.140
#Definition der Funktionen für Ausgleichsrechnung und Differntialquotienten
def f(x, a, b, c):
    return a * (x**2) + b * x + c

def f3(x, a, b, c, d):
    return a * (x**3) + b * (x**2) + c * x + d

def fstrich(x, a, b):
    return 2 * a * x + b

#Ausgleichsrechnung für Temperaturverläufe und speichern der Koeefizienten A, B, und C (5b) )
params1, covariance1 = curve_fit(f, ts, T1)
errors1 = np.sqrt(np.diag(covariance1))
a1 = ufloat(params1[0], errors1[0])
b1 = ufloat(params1[1], errors1[1])
c1 = ufloat(params1[2], errors1[2])

params2, covariance2 = curve_fit(f, ts, T2)
errors2 = np.sqrt(np.diag(covariance2))
a2 = ufloat(params2[0], errors2[0])
b2 = ufloat(params2[1], errors2[1])
c2 = ufloat(params2[2], errors2[2])

#Ausgleichsrechnung 3. Grades:
Params1, Covariance1 = curve_fit(f3, ts, T1)
Errors1 = np.sqrt(np.diag(Covariance1))
A1 = ufloat(Params1[0], Errors1[0])
B1 = ufloat(Params1[1], Errors1[1])
C1 = ufloat(Params1[2], Errors1[2])
D1 = ufloat(Params1[3], Errors1[3])

Params2, Covariance2 = curve_fit(f3, ts, T2)
Errors2 = np.sqrt(np.diag(Covariance2))
A2 = ufloat(Params2[0], Errors2[0])
B2 = ufloat(Params2[1], Errors2[1])
C2 = ufloat(Params2[2], Errors2[2])
D2 = ufloat(Params2[3], Errors2[3])

print("Aufgabe 5b), Parameter für Ax^3+Bx^2+Cx+D:")
print("T1")
print('A =', Params1[0], '±', Errors1[0])
print('B =', Params1[1], '±', Errors1[1])
print('C =', Params1[2], '±', Errors1[2])
print('C =', Params1[2], '±', Errors1[2])
print("T2")
print('A =', Params2[0], '±', Errors2[0])
print('B =', Params2[1], '±', Errors2[1])
print('C =', Params2[2], '±', Errors2[2])
print('D =', Params2[3], '±', Errors2[3])
print("")
#Berechnung der Quotienten (5c) )
Quotienten1 = f(zeiten, A1, B1, C1)
Quotienten2 = f(zeiten, A2, B2, C2)
print("Aufgabe 5c), Quotienten := Werte der Ableitung an den Stellen t")
print("Zeiten: ", zeiten)
print("Quotienten für T1: ", Quotienten1)
print("Quotienten für T2: ", Quotienten2)
print("")

#Berechnung der Güteziffern (5d) )
Leistung = ufloat(np.mean(N), stats.sem(N))
Temperatur1 = f3(zeiten, A1, B1, C1, D1)
#Versuch zu runden, klappt nicht...T1 = unp.around(f3(zeiten, A1, B1, C1, D1), decimals=3)
Temperatur2 = f3(zeiten, A2, B2, C2, D2)
nuemp = (cmapperat + cmwasser) * Quotienten1 * (1/Leistung)
nuid = Temperatur1 / (Temperatur1 - Temperatur2)
print("Aufgabe 5d), Güteziffern:")
print("Mittel der Leistung: ", Leistung)
print("Zeiten: ", zeiten)
print("T1 für die zeiten: ", Temperatur1)
print("T2 für die zeiten: ", Temperatur2)
print("Real: ", nuemp)
print("Ideal: ", nuid)
print("Mittelwert von N: ", Leistung)
print("")

#Ausgleichsrechnung zur Bestimmung der Verdampfungswärme L (5d) )
R = 8.314
M = 120.9 * (10**(-3))
def DDK(x, L, K):
    return L * x + K
Lparams, Lcovariance = curve_fit(DDK, 1/T1, np.log(Pb))
Lerrors = np.sqrt(np.diag(Lcovariance))
g = ufloat(Lparams[0], Lerrors[0])
h = ufloat(Lparams[1], Lerrors[1])
L = (- R * g) / M

#Bestimmung des Massendurchsatzes:
mdurchs = (cmapperat + cmwasser) * Quotienten2 * (1/L)
print("Aufgabe 5e, Verdampfungswärme und Massendurchsatz")
print("Verdampfungswärme: ", L)
print("Steigung: ", g)
print("Y-Achsenabschnitt: ", h)
print("Zeiten: ", zeiten)
print("Massendurchsatz: ", mdurchs)
print("")

#Bestimmung der Mechanischen Leistung, Formel für rho von Kevin/Mike
rho = (273.15*BDruck*rhogasnorm)/(100000*Temperatur1)
Nmech = 1/(kappa-1) * (ADruck * (BDruck/ADruck)**(1/kappa) - BDruck) * 1/rho * mdurchs
print("Zeiten: ", zeiten)
print("Dichte:", rho)
print("Mechanische Kompressorleistung: ", Nmech)
print("Danke, dass Sie sich für diese Auswertung entschieden haben.")

#print(h)

#Ausgaben


#Plot
x_plot = np.linspace(0, 18*60)

plt.figure(0)
plt.title("Temperaturverlauf in den Gefäßen mit Regression 2. Grades")

plt.ylabel('Temperatur (K)')
plt.xlabel("Zeit (s)")
plt.plot(ts, T1, 'r+', label="T1")
plt.plot(ts, T2, 'y+', label="T2")
plt.plot(x_plot, f(x_plot, *params1), 'b-', label='Regression für T1')
plt.plot(x_plot, f(x_plot, *params2), 'g-', label='Regression für T2')
plt.legend(loc="best")
plt.tight_layout
plt.savefig('Grad2.pdf')

plt.figure(1)
plt.title("Temperaturverlauf in den Gefäßen mit Regression 3. Grades")

plt.ylabel('Temperatur (K)')
plt.xlabel("Zeit (s)")
plt.plot(ts, T1, 'r+', label="T1")
plt.plot(ts, T2, 'y+', label="T2")
plt.plot(x_plot, f3(x_plot, *Params1), 'b-', label='Regression für T1')
plt.plot(x_plot, f3(x_plot, *Params2), 'g-', label='Regression für T2')
plt.legend(loc="best")
plt.tight_layout
plt.savefig('Grad3.pdf')

plt.figure(2)
plt.title("Temperaturverlauf in den Gefäßen über 18 Minuten gemessen")

plt.ylabel('Temperatur (K)')
plt.xlabel("Zeit (s)")
plt.plot(ts, T1, 'r+', label="T1")
plt.plot(ts, T2, 'b+', label="T2")
plt.legend(loc="best")
plt.tight_layout
plt.savefig('verlauf.pdf')
