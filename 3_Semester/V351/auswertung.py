import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 15

#Einlesen der Daten
f_D_R_in_kHz, U_Rechteck, U_Dreieck = np.genfromtxt('rechteck_dreieck.txt', unpack=True)
f_S_in_kHz, U_Saegezahn = np.genfromtxt('saegezahn.txt', unpack=True)

Vpp = 10
R = 50
Frequenz_Funktionsgenerator = 10000

f_D_R = f_D_R_in_kHz * 1000
f_S = f_S_in_kHz * 1000

U_R_erw = 8.88 * 10/f_D_R_in_kHz
U_D_erw = 5.6 * 1/(f_D_R_in_kHz/10)**2
U_S_erw = 4.4 * 10/f_S_in_kHz

np.savetxt('Tabelle_R_D.txt', np.column_stack([f_D_R_in_kHz, U_Rechteck, U_R_erw, U_Dreieck, U_D_erw]), fmt="%.2f")
np.savetxt('Tabelle_S.txt', np.column_stack([f_S_in_kHz, U_Saegezahn, U_S_erw]), fmt="%.2f" )
