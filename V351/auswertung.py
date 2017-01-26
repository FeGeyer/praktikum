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
f_S_in_kHz = np.genfromtxt('saegezahn.txt', unpack=True)

Vpp = 10
R = 50
Frequenz_Funktionsgenerator = 10000

f_D_R = f_D_R_in_kHz * 1000
f_S = f_S_in_kHz * 1000
