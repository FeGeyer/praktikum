import numpy as np
import matplotlib.pyplot as plt
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

print("G1: ", g1)
print("g2: ", g2)
g = g1/g2
print("Verhältnis 1/2: ", g)

print('Horizontalkomponente 1: ', bB1)
print('Horizontalkomponente 2: ', bB2)

x = np.linspace(0, 10**6, 100000)
plt.plot(RFfrequenz, B1, 'r.', label="Erstes Isotop")
plt.plot(x, linear(x, *paramsB1), 'k--', label="Fit 1")
plt.plot(RFfrequenz, B2, 'b.', label="Zweites Isotop")
plt.plot(x, linear(x, *paramsB2), 'g--', label="Fit 2")
plt.xlabel(r'RF-Frequenz / Hz')
plt.ylabel(r'Horizontalkomponente von B / T')
plt.tight_layout()
plt.legend(loc='best')
plt.savefig('Auswertung/Plots/Bfeld.pdf')


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

print("b Parameter1: ", bPeriode1)

# fit for second period
paramsPeriode2, covPeriode2 = curve_fit(hyperbel, amplitude[1:], periode2[1:])
errorsPeriode2 = np.sqrt(np.diag(covPeriode2))
aPeriode2 = ufloat(paramsPeriode2[0], errorsPeriode2[0])
bPeriode2 = ufloat(paramsPeriode2[1], errorsPeriode2[1])
cPeriode2 = ufloat(paramsPeriode2[2], errorsPeriode2[2])

print('b Parameter2: ', bPeriode2)
print("Verhältnis: ", bPeriode2/bPeriode1)

# plot this nice plots
x = np.linspace(1, 10.5, 10000)
plt.plot(amplitude, periode1, 'b.', label="Periode 1")
plt.plot(amplitude[1:], periode2[1:], 'r.', label="Periode 2")
plt.plot(x, hyperbel(x, *paramsPeriode1), 'k--',label="Fit Periode 1")
plt.plot(x, hyperbel(x, *paramsPeriode2), 'g--',label="Fit Periode 2")
plt.xlabel('RF-Amplitude / V')
plt.ylabel(r'Periodendauer/ ms')
plt.tight_layout()
plt.legend(loc='best')
plt.savefig('Auswertung/Plots/Perioden.pdf')
