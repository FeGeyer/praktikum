import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

# Impuls-Echo-Verfahren
# Einlesen der Daten
IlaengeMM, IampliV, ItMikros, IdeltatMikros, ITGCdB = np.genfromtxt('impuls.txt', unpack=True)
DlaengeMM, DampliV, DtMikros, DTGCdB = np.genfromtxt('durchschall.txt', unpack=True)
Apeak, AtMikros, AdeltatMikros = np.genfromtxt('auge.txt', unpack=True)
Cpeak, CtMikros = np.genfromtxt('cepstrum.txt', unpack=True)
print(AtMikros)
print(IampliV)
