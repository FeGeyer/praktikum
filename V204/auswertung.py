import numpy as np

Nr, s_statisch, T1_statisch, T4_statisch, T5_statisch, T8_statisch = np.genfromtxt('statisch.txt', unpack=True)

print(s_statisch)
print(T1_statisch)
