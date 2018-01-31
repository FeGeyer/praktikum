import numpy as np
from uncertainties import ufloat
from scipy.stats import stats

D_theta, R1, R2, R3, R4, R5 = np.genfromtxt('glas.txt', unpack=True)

T = 1 * 10**-3
lam = 633 * 10**(-9)
theta_0 = 138 * np.pi / 180
theta = D_theta * (np.pi)/180


def n(lam, N, T, t):
    return (1-(lam * N)/(2*T*0.175*t))**(-1)


n1 = ufloat(np.mean(n(lam, R1, T, theta)), stats.sem(n(lam, R1, T, theta)))
n2 = ufloat(np.mean(n(lam, R2, T, theta)), stats.sem(n(lam, R2, T, theta)))
n3 = ufloat(np.mean(n(lam, R3, T, theta)), stats.sem(n(lam, R3, T, theta)))
n4 = ufloat(np.mean(n(lam, R4, T, theta)), stats.sem(n(lam, R4, T, theta)))
n5 = ufloat(np.mean(n(lam, R5, T, theta)), stats.sem(n(lam, R5, T, theta)))


n_array = [n1, n2, n3, n4, n5]
n_best = np.mean(n_array)

print(n_best)
