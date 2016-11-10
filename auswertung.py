#import matplotlib.pyplot as plt
#plt.rcParams['figure.figsize'] = (10, 8)
#plt.rcParams['font.size'] = 16

import numpy as np
from scipy.optimize import curve_fit
from scipy import stats
import uncertainties.unumpy as unp
from uncertainties import ufloat


T_1r, T_2r, T_pr, T_mr, Tr = np.genfromtxt('data_60cm', unpack=True)
T_1 = ufloat(np.mean(T_1r/5), stats.sem(T_1r/5))
T_2 = ufloat(np.mean(T_2r/5), stats.sem(T_2r/5))
T_p = ufloat(np.mean(T_pr/5), stats.sem(T_pr/5))
T_m = ufloat(np.mean(T_mr/5), stats.sem(T_mr/5))
T = ufloat(np.mean(Tr), stats.sem(Tr))
print(T_1)
print(T_2)
print(T_p)
print(T_m)
print(T)
