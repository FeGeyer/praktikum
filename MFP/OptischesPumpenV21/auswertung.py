import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
import pandas as pd
from scipy.optimize import curve_fit
from uncertainties import ufloat
from scipy import constants

# Use latex fonts and text
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Define required functions
def linear(x, m, b):
    return x*m + b

def hyperbel(x, a, b, c):
    return a + b/(x-c)

def Bfeld(I, N, R):
    return constants.mu_0 * (8*I*N)/(np.sqrt(125)*R)

def exponential(x, a, b, c):
    return -a**(x*b) - c

bohr = constants.value('Bohr magneton')
# load data
RFfrequenz, Sweep1, Horizontal1, Sweep2, Horizontal2 = np.genfromtxt('Auswertung/Daten/daten_sweep.txt', unpack=True)


RFfrequenz = RFfrequenz*10**(3)
# scale in ampere
I1_Sweep = Sweep1*0.1
I2_Sweep = Sweep2*0.1

I1_Horizontal = Horizontal1*0.3
I2_Horizontal = Horizontal2*0.3

# Constants
N_Sweep = 11
R_Sweep = 16.39*10**(-2)

N_Horizontal = 154
R_Horizontal = 15.79*10**(-2)

N_Vertical = 20
R_Vertical = 11.735*10**(-2)

# calculate magnetic field
B1_Sweep = Bfeld(I1_Sweep, N_Sweep, R_Sweep)
B1N_Horizontal = Bfeld(I1_Horizontal, N_Horizontal, R_Horizontal)
B1 = (B1_Sweep + B1N_Horizontal)

B2_Sweep = Bfeld(I2_Sweep, N_Sweep, R_Sweep)
B2_Horizontal = Bfeld(I2_Horizontal, N_Horizontal, R_Horizontal)
B2 = (B2_Sweep + B2_Horizontal)

# fit magnetic field against frequency
paramsB1, covB1 = curve_fit(linear, RFfrequenz, B1)
errorsB1 = np.sqrt(np.diag(covB1))
mB1 = ufloat(paramsB1[0], errorsB1[0])
bB1 = ufloat(paramsB1[1], errorsB1[1])

paramsB2, covB2 = curve_fit(linear, RFfrequenz, B2)
errorsB2 = np.sqrt(np.diag(covB2))
mB2 = ufloat(paramsB2[0], errorsB2[0])
bB2 = ufloat(paramsB2[1], errorsB2[1])

# calculate g_F
g1 = constants.h/(mB1 * bohr)
g2 = constants.h/(mB2 * bohr)

print('Landésche gF')
print("g1: ", g1)
print("g2: ", g2)
g = g1/g2
print("Verhältnis 1/2: ", g)
print('')

print('---------------------------------------------------------------------')
print('Horizontalkomponenten')
print('Steigung 1: ', mB1)
print('Steigung 2: ', mB2)
print('Horizontalkomponente 1: ', bB1)
print('Horizontalkomponente 2: ', bB2)
print('')

x = np.linspace(0, 10**6, 100000)
plt.plot(RFfrequenz, B1, 'r.', markersize=10, fillstyle='none', label="Erstes Isotop")
plt.plot(x, linear(x, *paramsB1), 'r-', label="Fit 1")
plt.plot(RFfrequenz, B2, 'b.', markersize=10, fillstyle='none', label="Zweites Isotop")
plt.plot(x, linear(x, *paramsB2), 'b--', label="Fit 2")
plt.xticks([0, 200000, 400000, 600000, 800000, 1000000], ['0', '0.2', '0.4', '0.6', '0.8', '1'])
plt.yticks([0, 0.00005, 0.00010, 0.00015, 0.00020], ['0', '50', '100', '150', '200'])
plt.xlabel(r'RF-Frequenz / MHz')
plt.ylabel(r'Horizontalkomponente von B / $\mu$T')
plt.tight_layout()
plt.legend(loc='best')
plt.savefig('Auswertung/Plots/Bfeld.pdf')
plt.clf()


amplitude, periode1, periode2 = np.genfromtxt('Auswertung/Daten/daten_periode.txt', unpack=True)

# scale to ms
periode1 = periode1*10**(-3)
periode2 = periode2*10**(-3)

# fit for first period
paramsPeriode1, covPeriode1 = curve_fit(hyperbel, amplitude[1:], periode1[1:])
errorsPeriode1 = np.sqrt(np.diag(covPeriode1))
aPeriode1 = ufloat(paramsPeriode1[0], errorsPeriode1[0])
bPeriode1 = ufloat(paramsPeriode1[1], errorsPeriode1[1])
cPeriode1 = ufloat(paramsPeriode1[2], errorsPeriode1[2])

print('--------------------------------------------------------------')
print('Hyperbel b Parameter')
print('a und c 1: ', aPeriode1, cPeriode1)
print("b Parameter1: ", bPeriode1)

# fit for second period
paramsPeriode2, covPeriode2 = curve_fit(hyperbel, amplitude[1:], periode2[1:])
errorsPeriode2 = np.sqrt(np.diag(covPeriode2))
aPeriode2 = ufloat(paramsPeriode2[0], errorsPeriode2[0])
bPeriode2 = ufloat(paramsPeriode2[1], errorsPeriode2[1])
cPeriode2 = ufloat(paramsPeriode2[2], errorsPeriode2[2])

print('a und c 2: ', aPeriode2, cPeriode2)
print('b Parameter2: ', bPeriode2)
print("Verhältnis 1/2: ", bPeriode1/bPeriode2)
print("Verhältnis 2/1: ", bPeriode2/bPeriode1)
print('')

# plot this nice plots
x = np.linspace(1, 10.5, 10000)
plt.plot(amplitude, periode1, 'b.', markersize=10, fillstyle='none', label="Periode 1")
plt.plot(amplitude[1:], periode2[1:], 'r.', markersize=10, fillstyle='none', label="Periode 2")
plt.plot(x, hyperbel(x, *paramsPeriode1), 'b--',label="Fit Periode 1")
plt.plot(x, hyperbel(x, *paramsPeriode2), 'r-',label="Fit Periode 2")
plt.xlabel('RF-Amplitude / V')
plt.ylabel(r'Periodendauer/ ms')
plt.tight_layout()
plt.legend(loc='best')
plt.savefig('Auswertung/Plots/Perioden.pdf')
plt.clf()

# calculate nuclear spin
J = 0.5
S = 0.5
L = 0
gJ = (3.0023 * (J**2 + J) + 1.0023 * ((S**2 + S)
      - (L**2 + L))) / (2 * (J**2 + J))
I1 = gJ / (4 * g1) - 1 + unp.sqrt((gJ / (4 * g1) - 1)**2
                                  + 3 * gJ / (4 * g1) - 3 / 4)
I2 = gJ / (4 * g2) - 1 + unp.sqrt((gJ / (4 * g2) - 1)**2
+ 3 * gJ / (4 * g2) - 3 / 4)

print('--------------------------------------------------------------')
print('Kernspins')
print('Kernspin 1: ', I1)
print('Kernspin 2:', I2)
print('')

# assessment quadratic zeeman
U1 = g1*bohr*np.max(B1)+g1**2*bohr**2*np.max(B1)**2*(1-2*2)/(4.53e-24)
U2 = g2*bohr*np.max(B2)+g2**2*bohr**2*np.max(B2)**2*(1-2*3)/(2.01e-24)

print('------------------------------------------------------------')
print('Quadratische Zeeman-Aufspaltung')
print('Maximales BFeld1: ', np.round(np.max(B1)*10**6, 2))
print('Maximales BFeld2: ', np.round(np.max(B2)*10**6, 2))
print('Quadratische Zeeman-Aufspaltung 1 in eV: ', U1/constants.e)
print('Quadratische Zeeman-Aufspaltung 2 in eV: ', U2/constants.e)

# Exponential fit
#exp1 = pd.DataFrame()
#exp2 = pd.DataFrame()
#exp1 = pd.read_csv('Auswertung/Daten/1.CSV')
#exp2 = pd.read_csv('Auswertung/Daten/TEK0000.CSV')
#
#print(exp1)

#x1, y1 = np.genfromtxt('Auswertung/Daten/1.txt', unpack=True)
#x2, y2 = np.genfromtxt('Auswertung/Daten/2.txt', unpack=True)
#x1 = x1[979:1582]
#y1 = y1[979:1582]
#
#paramsExp1, covExp1 = curve_fit(exponential, x1, y1)
#
#x = np.linspace(-0.0002, 0.075, 10000)
#plt.plot(x1, y1, label="Daten")
#plt.plot(x, exponential(x, *paramsExp1), 'r--', label='Fit')
#plt.legend(loc="best")
##plt.xscale('log')
##plt.plot(x2, y2)
#plt.savefig('Auswertung/Plots/exp1.pdf')
