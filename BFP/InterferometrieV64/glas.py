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


n1_mean = ufloat(np.mean(n(lam, R1, T, theta)), stats.sem(n(lam, R1, T, theta)))
n2_mean = ufloat(np.mean(n(lam, R2, T, theta)), stats.sem(n(lam, R2, T, theta)))
n3_mean = ufloat(np.mean(n(lam, R3, T, theta)), stats.sem(n(lam, R3, T, theta)))
n4_mean = ufloat(np.mean(n(lam, R4, T, theta)), stats.sem(n(lam, R4, T, theta)))
n5_mean = ufloat(np.mean(n(lam, R5, T, theta)), stats.sem(n(lam, R5, T, theta)))

n1 = n(lam, R1, T, theta)
n2 = n(lam, R2, T, theta)
n3 = n(lam, R3, T, theta)
n4 = n(lam, R4, T, theta)
n5 = n(lam, R5, T, theta)

# np.append(n1, n1_mean)
# np.append(n2, n2_mean)
# np.append(n3, n3_mean)
# np.append(n4, n4_mean)
# np.append(n5, n5_mean)
# np.append(theta, np.nan)

print(n1)

np.savetxt('TexTabellen/glas.txt', np.column_stack([
        theta * 180/np.pi,
        n1,
        n2,
        n3,
        n4,
        n5
        ]), delimiter=' & ', newline=r' \\'+'\n',
        fmt='%.0f & %.2f & %.2f & %.2f & %.2f & %.2f')


n_array = [n1_mean, n2_mean, n3_mean, n4_mean, n5_mean]
n_best = np.mean(n_array)

print(n_best)
