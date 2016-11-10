#import matplotlib.pyplot as plt
#plt.rcParams['figure.figsize'] = (10, 8)
#plt.rcParams['font.size'] = 16

import numpy as np
from scipy.optimize import curve_fit
from scipy import stats
import uncertainties.unumpy as unp
from uncertainties import ufloat


AT_1r, AT_2r, AT_pr, AT_mr, ATr = np.genfromtxt('data_60cm', unpack=True)
AT_sr = np.genfromtxt('data_60_Ts', unpack=True)
AT_1 = ufloat(np.mean(AT_1r/5), stats.sem(AT_1r/5))
AT_2 = ufloat(np.mean(AT_2r/5), stats.sem(AT_2r/5))
AT_p = ufloat(np.mean(AT_pr/5), stats.sem(AT_pr/5))
AT_m = ufloat(np.mean(AT_mr/5), stats.sem(AT_mr/5))
AT = ufloat(np.mean(ATr), stats.sem(ATr))
AT_s = ufloat(np.mean(AT_sr/5), stats.sem(AT_sr/5))

BT_1r, BT_2r, BT_pr, BT_mr, BTr = np.genfromtxt('data_75cm', unpack=True)
BT_sr = np.genfromtxt('data_75_Ts', unpack=True)
BT_1 = ufloat(np.mean(BT_1r/5), stats.sem(BT_1r/5))
BT_2 = ufloat(np.mean(BT_2r/5), stats.sem(BT_2r/5))
BT_p = ufloat(np.mean(BT_pr/5), stats.sem(BT_pr/5))
BT_m = ufloat(np.mean(BT_mr/5), stats.sem(BT_mr/5))
BT = ufloat(np.mean(BTr), stats.sem(BTr))
BT_s = ufloat(np.mean(BT_sr/5), stats.sem(BT_sr/5))

print("1T_1:", AT_1)
print("1T_2:", AT_2)
print("1T_p:", AT_p)
print("1T_m:", AT_m)
print("1T:", AT)
print("1T_s:", AT_s)

print("2T_1:", BT_1)
print("2T_2:", BT_2)
print("2T_p:", BT_p)
print("2T_m:", BT_m)
print("2T:", BT)
print("2T_s:", BT_s)
