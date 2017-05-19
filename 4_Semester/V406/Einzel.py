import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = (14, 8)
plt.rcParams['font.size'] = 16

# Messwerte: x (Pos. der Photodiode), y (Intensität des Lasers)
x, y = np.genfromtxt('einzel.txt', unpack=True)

b_mat = 0.075           # Materialwert der Spaltbreite
y *= 1               # Umrechnung der Intensität in nA
y -= 0.155             # Abzug des Dunkelstroms
L = 1050               # Abstand zwischen Spalt und Detektor
l = 6.33*(10**(-4))     # Wellenlänge des Lasers in mm


# Fit der Theoriekurve an die Messwerte
def f(x, b, A, B):
    return (A**2)*(b**2)*(l/(np.pi*b*np.sin((x-B)/L)))**2*(np.sin((np.pi*b*np.sin((x-B)/L))/l))**2

#                   geschaetzte Parameter: Spaltbr, max Ampl.  x-Abst. 0-Max
params, covariance = curve_fit(f, x, y, p0=[0.075, np.max(y), 30.1])

errors = np.sqrt(np.diag(covariance))

print('b =', params[0], '±', errors[0])
print('A =', params[1], '±', errors[1])
print('B =', params[2], '±', errors[2])
#print('L =', params[3], '±', errors[3])

x_plot = np.linspace(9, 52, 500)

plt.plot(x, y, 'rx', label="Messwerte")
#plt.plot(x_plot, (A**2)*(b**2)*(l/(np.pi*b*np.sin((x_plot-B)/L)))**2*(np.sin((np.pi*b*np.sin((x_plot-B)/L))/l))**2, 'r-', label='theorie')
plt.plot(x_plot, f(x_plot, *params), 'b-', label='nicht-linearer Fit')

plt.xlim(min(x)-2, max(x)+2)
plt.ylim(min(y)-20, max(y)+50)

plt.xlabel(r'd \ mm')
plt.ylabel(r'I \ nA')

plt.legend(loc="best")
plt.savefig('einzel_test.png')

# Abweichung des berechneten Werts von dem Materialwert
print(abs(100 - (b_mat/params[0] * 100)))
