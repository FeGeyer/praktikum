import numpy as np
import matplotlib.pyplot as plt

# Einlesen der Daten
# Leermessung
daten = np.genfromtxt('../Daten/Wuerfel%1.0u/I%1.0u.Spe' % (1, 1), unpack=True)

# for i in range(3):
#     for j in range(12):
#         daten = np.genfromtxt('Daten/Wuerfel%1.0u/I%1.0u' % (i+1, i+1))
#
#daten_1_1 = np.genfromtxt('Daten/Wuerfel1/I1.Spe')
#daten_1_1 = np.genfromtxt('Daten/Wuerfel1/I2.Spe')
#daten_1_1 = np.genfromtxt('Daten/Wuerfel1/I3.Spe')
#daten_1_1 = np.genfromtxt('Daten/Wuerfel1/I4.Spe')
#daten_1_1 = np.genfromtxt('Daten/Wuerfel1/I5.Spe')
#daten_1_1 = np.genfromtxt('Daten/Wuerfel1/I6.Spe')
#daten_1_1 = np.genfromtxt('Daten/Wuerfel1/I7.Spe')
#daten_1_1 = np.genfromtxt('Daten/Wuerfel1/I8.Spe')
#daten_1_1 = np.genfromtxt('Daten/Wuerfel1/I9.Spe')
#daten_1_1 = np.genfromtxt('Daten/Wuerfel1/I1.Spe')
#daten_1_1 = np.genfromtxt('Daten/Wuerfel1/I1.Spe')
#daten_1_1 = np.genfromtxt('Daten/Wuerfel1/I1.Spe')
i = 1
print('Wuerfel %1.0u' % (2))
