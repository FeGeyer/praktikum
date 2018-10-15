import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

# Einlesen der Daten
# Umrechnungsfaktor für Energie
umrechnung = 9.323944
# Leermessung

daten = np.genfromtxt('../Daten/Quellspektrum/Quellspektrum_300.Spe', unpack=True)
x = np.linspace(0, 511, 512)
plt.plot(x*umrechnung, daten/3, 'r--',label="Daten")
plt.axvline(662, color='black',label=r"$662 \ \mathrm{keV}$")
plt.xlim(0, 1000)
plt.ylim(-50, 3500)
plt.yticks([0, 500, 1000, 1500, 2000, 2500, 3000], ['0', '0,5', '1', '1,5', '2', '2,5', '3'])
plt.xlabel(r"$\mathrm{Energie} \ / \ \mathrm{keV}$")
plt.ylabel(r'$\mathrm{Anzahl \ 1000 \ Ereignisse} \ / \ 100 \, \mathrm{s}$')
plt.legend(loc='best')
plt.savefig('quellmessung.pdf')
plt.clf()
# Würfel 1

for i in range(12):
    daten = np.genfromtxt('../Daten/Wuerfel%1.0u/I%1.0u.Spe' % (1, i+1), unpack=True)
    x = np.linspace(0, 511, 512)
    plt.plot(x*umrechnung, daten, 'r--',label="Daten")
    plt.axvline(662, color='black',label=r"$662 \ \mathrm{keV}$")
    plt.xlim(0, 1000)
    plt.ylim(-50, 3500)
    plt.yticks([0, 500, 1000, 1500, 2000, 2500, 3000], ['0', '0,5', '1', '1,5', '2', '2,5', '3'])
    plt.xlabel(r"$\mathrm{Energie} \ / \ \mathrm{keV}$")
    plt.ylabel(r'$\mathrm{Anzahl \ 1000 \ Ereignisse} \ / \ 100 \, \mathrm{s}$')
    plt.legend(loc='best')
    plt.savefig('wuerfel1_%1.0u.pdf' % (i+1))
    plt.clf()
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
print('Alles ausgeführt!')
