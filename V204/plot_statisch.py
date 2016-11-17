import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 16


Nr, s_statisch, T1_statisch, T4_statisch, T5_statisch, T8_statisch = np.genfromtxt('Statisch.txt', unpack=True)


plt.ylabel("Temperatur in Celsius")
plt.xlabel("Zeit in s")
plt.plot(s_statisch, T1_statisch, 'r--', label="Messing breit")
plt.plot(s_statisch, T4_statisch, 'y--', label="Messing dünn")
plt.plot(s_statisch, T5_statisch, 'g--', label="Aluminium")
plt.plot(s_statisch, T8_statisch, 'b--', label="Edelstahl")
#plt.plot(x_plot, f(x_plot, *params), 'b-', label='Regression', linewidth=3)
plt.legend(loc="best")
plt.savefig('Statisch.pdf')
