#import matplotlib.pyplot as plt
#plt.rcParams['figure.figsize'] = (10, 8)
#plt.rcParams['font.size'] = 16

import numpy as np
from scipy.optimize import curve_fit
from scipy import stats
import scipy.constants as const
import uncertainties.unumpy as unp
from uncertainties import ufloat

x=1
AT_1r, AT_2r, AT_pr, AT_mr, ATr = np.genfromtxt('data_60cm', unpack=True)
AT_sr = np.genfromtxt('data_60_Ts', unpack=True)
AT_1 = ufloat(np.mean(AT_1r/5), stats.sem(AT_1r/5))
AT_2 = ufloat(np.mean(AT_2r/5), stats.sem(AT_2r/5))
AT_p = ufloat(np.mean(AT_pr/5), stats.sem(AT_pr/5))
AT_m = ufloat(np.mean(AT_mr/5), stats.sem(AT_mr/5))
AT = ufloat(np.mean(ATr), stats.sem(ATr))
AT_s = ufloat(np.mean(AT_sr/5), stats.sem(AT_sr/5))
Al = ufloat(0.6, 0.003)
X = ufloat(x, 7)
print(X)

BT_1r, BT_2r, BT_pr, BT_mr, BTr = np.genfromtxt('data_75cm', unpack=True)
BT_sr = np.genfromtxt('data_75_Ts', unpack=True)
BT_1 = ufloat(np.mean(BT_1r/5), stats.sem(BT_1r/5))
BT_2 = ufloat(np.mean(BT_2r/5), stats.sem(BT_2r/5))
BT_p = ufloat(np.mean(BT_pr/5), stats.sem(BT_pr/5))
BT_m = ufloat(np.mean(BT_mr/5), stats.sem(BT_mr/5))
BT = ufloat(np.mean(BTr), stats.sem(BTr))
BT_s = ufloat(np.mean(BT_sr/5), stats.sem(BT_sr/5))
Bl = ufloat(0.75, 0.003)

Aomegaplus = ((2 * np.pi) /  AT_p)
Aomegaminus = (2 * np.pi) / AT_m
Aomegas = (2 * np.pi) / AT_s
AK = ((Aomegaminus ** 2) - (Aomegaplus ** 2)) / ((Aomegaminus ** 2) + (Aomegaplus ** 2))

Bomegaplus = (2 * np.pi) /  BT_p
Bomegaminus = (2 * np.pi) / BT_m
Bomegas = (2 * np.pi) / BT_s
BK = ((Bomegaminus ** 2) - (Bomegaplus ** 2)) / ((Bomegaminus ** 2) + (Bomegaplus ** 2))

RechnerischSchw1 = (AT_p*AT_m)/(AT_p-AT_m)
RechnerischSchw2 = (BT_p*BT_m)/(BT_p-BT_m)

a=np.around(AT_pr/5, decimals=2)
b=np.around(AT_mr/5, decimals=2)
x=np.around(BT_pr/5, decimals=2)
y=np.around(BT_mr/5, decimals=2)
print(a)
# np.savetxt('tabelle.tex', np.column_stack([a, b]), fmt="%.2f")
# np.savetxt('tabelle2.tex', np.column_stack([x, y]), fmt="%.2f")
np.savetxt('1.txt', np.column_stack([a, b, x, y]), fmt="%.2f")
print('Länge: 60cm')
print('T_1: ', AT_p)
print('T_2: ', AT_m)
print('Omega+: ', Aomegaplus)
print('Omega-: ', Aomegaminus)
print('Kopplung: ', AK)
print('Berechnete Schwebungsdauer: ', RechnerischSchw1)
print('Gemessene Einzelschwebungsdauer: ', AT)
print('Gemessene Mehrfachschwingungsdauer: ', AT_s, Aomegas)
print(' ')
print('Länge: 75cm')
print('T_1: ', BT_p)
print('T_2: ', BT_m)
print('Omega+: ', Bomegaplus)
print('Omega-: ', Bomegaminus)
print('Kopplung: ', BK)
print('Berechnete Schwebungsdauer: ', RechnerischSchw2)
print('Gemessene Einzelschwebungsdauer: ', BT)
print('Gemessene Mehrfachschwingungsdauer: ', BT_s, Bomegas)
