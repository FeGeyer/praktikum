import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 10


tmin, T1, Pb, T2, Pa, N = np.genfromtxt('data.txt', unpack=True)
ts = tmin*60
T1 = T1+273.15
T2 = T2+273.15

def f(x, a, b, c):
    return a * (x**2) + b * x + c
params1, covariance1 = curve_fit(f, ts, T1)
params2, covariance2 = curve_fit(f, ts, T2)
errors1 = np.sqrt(np.diag(covariance1))
errors2 = np.sqrt(np.diag(covariance2))
print("T1")
print('a =', params1[0], '±', errors1[0])
print('b =', params1[1], '±', errors1[1])
print('c =', params1[2], '±', errors1[2])
print("T2")
print('a =', params2[0], '±', errors2[0])
print('b =', params2[1], '±', errors2[1])
print('c =', params2[2], '±', errors2[2])

plt.title("Temperaturverlauf in den Gefäßen")

x_plot = np.linspace(0, 18*60)

plt.ylabel('Temperatur (K)')
plt.xlabel("Zeit (s)")
plt.plot(ts, T1, 'r+', label="T1")
plt.plot(ts, T2, 'y+', label="T2")
plt.plot(x_plot, f(x_plot, *params1), 'b-', label='Regression', linewidth=3)
plt.plot(x_plot, f(x_plot, *params2), 'b-', label='Regression', linewidth=3)
plt.legend(loc="best")
plt.savefig('Tverl.pdf')

plt.tight_layout()
