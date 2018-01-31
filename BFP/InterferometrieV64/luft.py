import numpy as np
from uncertainties import ufloat
from uncertainties import unumpy
from scipy.stats import stats

N = np.genfromtxt('luft.txt', unpack=True)

N = unumpy.uarray(N, 4)

def n(N, lam, L):
    return N * lam / L + 1


lam = 633 * 10**(-9)
L = 0.1

n_array = n(N, lam, L)

n_best = np.mean(n_array)
print(n_best)
