import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
import scipy.constants as const
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

deltad, z_well = np.genfromtxt('wellenlaenge.txt', unpack=True)
deltaP, z_Luft, z_CO2 = np.genfromtxt('brechnungsind.txt', unpack=True)

uebers = 5.017
b = 50*10**(-3)
p0 = 1.0132
T0 = 273.15
T = 293.15
l = 665 * 10 ** (-9)

deltaP = p0 - deltaP

deltad = deltad * 10**(-3)/uebers

wl = 2*deltad/z_well
lambdaLaser = ufloat(np.mean(wl), stats.sem(wl))
print("Lambda:",lambdaLaser*10**9)

deltanLuft = z_Luft*l/(2*b)
deltanCO2 = z_CO2*l/(2*b)

pstrich = p0 - deltaP

nLuft = 1 + deltanLuft * (T/T0) * (p0/pstrich)
nCO2 = 1 + deltanCO2 * (T/T0) * (p0/pstrich)

nCO2Mean = ufloat(np.mean(unp.nominal_values(nCO2)), stats.sem(nCO2))
nLuftMean = ufloat(np.mean(unp.nominal_values(nLuft)), stats.sem(nLuft))
print(nLuftMean)
print(nCO2Mean)

np.savetxt('Wellenl.txt', np.column_stack([deltad * 10**3, z_well, wl*10**9]) , fmt="%.2f")
np.savetxt('Luft.txt', np.column_stack([deltaP, z_Luft, nLuft]) , fmt="%.6f")
np.savetxt('CO2.txt', np.column_stack([deltaP, z_CO2, nCO2]) , fmt="%.6f")
