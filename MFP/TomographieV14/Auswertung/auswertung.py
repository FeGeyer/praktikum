import numpy as np
import matplotlib.pyplot as plt
from uncertainties import unumpy
from uncertainties import ufloat
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18


#   Wegstrecke
d = 0.01    # m
# Umrechnungsfaktor für Energie
umrechnung = 9.323944

# Liste für Wüfelnummern
liste = [1, 2, 4]

# Messzeiten
t_100 = 100     # Würfel 1 und 2
t_300 = 300     # Würfel 4

# Array mit Peaks
peaks = np.array([])

# Array mit Counts
counts = np.array([])

# Leermessung

daten = np.genfromtxt('../Daten/Quellspektrum/Quellspektrum_echt_300.Spe',
                      unpack=True)
x = np.linspace(0, 511, 512)
plt.plot(x*umrechnung, daten/3, 'r--', label="Daten")
plt.axvline(662, color='black', label=r"$662 \ \mathrm{keV}$")
plt.xlim(0, 1000)
plt.ylim(-50, 3700)
plt.yticks([0, 500, 1000, 1500, 2000, 2500, 3000],
           ['0', '0,5', '1', '1,5', '2', '2,5', '3'])
plt.xlabel(r"$\mathrm{Energie} \ / \ \mathrm{keV}$")
plt.ylabel(r'$\mathrm{Anzahl \ 1000 \ Ereignisse} \ / \ 100 \, \mathrm{s}$')
plt.legend(loc='best')
plt.savefig('Plots/Quellmessung/quellmessung.pdf')
plt.clf()

# Plots für die Würfel erstellen

for j in liste:
    for i in range(12):
        daten = np.genfromtxt('../Daten/Wuerfel%1.0u/I%1.0u.Spe' % (j, i+1),
                              unpack=True)
        x = np.linspace(0, 511, 512)
        plt.plot(x*umrechnung, daten, 'r--', label="Daten")
        plt.axvline(662, color='black', label=r"$662 \ \mathrm{keV}$")
        plt.xlim(0, 1000)
        plt.ylim(-50, 3500)
        plt.yticks([0, 500, 1000, 1500, 2000, 2500, 3000],
                   ['0', '0,5', '1', '1,5', '2', '2,5', '3'])
        plt.xlabel(r"$\mathrm{Energie} \ / \ \mathrm{keV}$")
        plt.ylabel(r'$\mathrm{Anzahl \ 1000 \ Ereignisse} \ / \ 100 \, \mathrm{s}$')
        plt.legend(loc='best')
        plt.savefig('Plots/Wuerfel%1.0u/wuerfel%1.0u_%1.0u.pdf' % (j, j, i+1))
        plt.clf()
        peaks = np.append(peaks, 71)
        counts = np.append(counts, daten[71])


# Umrechnen der Peaks in Energie
peaks = peaks * umrechnung

peaks_1 = peaks[0:12]
peaks_2 = peaks[12:24]
peaks_4 = peaks[24:36]

counts_1 = counts[0:12]
counts_2 = counts[12:24]
counts_4 = counts[24:36]

# Berechnung von I_0
I_0_300 = 10762
t_i_0 = 300     # Sekunden

# Würfelmatrix
a = np.sqrt(2)
A = np.array([[0, a, 0, a, 0, 0, 0, 0, 0],  # I1
              [0, 0, a, 0, a, 0, a, 0, 0],
              [0, 0, 0, 0, 0, a, 0, a, 0],
              [1, 1, 1, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 1, 1, 1, 0, 0, 0],  # I5
              [0, 0, 0, 0, 0, 0, 1, 1, 1],
              [0, a, 0, 0, 0, a, 0, 0, 0],
              [a, 0, 0, 0, a, 0, 0, 0, a],
              [0, 0, 0, a, 0, 0, 0, a, 0],
              [0, 0, 1, 0, 0, 1, 0, 0, 1],  # I10
              [0, 1, 0, 0, 1, 0, 0, 1, 0],
              [1, 0, 0, 1, 0, 0, 1, 0, 0]])


def kleinsteQuadrate(y, W, A):
    temp = np.dot(np.linalg.inv(np.dot(A.T, np.dot(W, A))), A.T)
    a = np.dot(temp, np.dot(W, y))
    a_err = np.linalg.inv(np.dot(A.T, np.dot(W, A)))
    return a, np.sqrt(np.diag(a_err))


# Raten berechnen
# Nullrate
I_0 = I_0_300/t_i_0
err_I_0 = np.sqrt(I_0)/t_i_0
print("Nullrate: ", I_0)
print("Fehler auf der Nullrate: ", np.round(err_I_0, 3))

# Würfel 1
rate_1 = counts_1/t_100
err_rate_1 = np.sqrt(rate_1)/t_100

# Würfel 2
rate_2 = counts_2/t_100
err_rate_2 = np.sqrt(rate_2)/t_100

# Würfel 4
rate_4 = counts_4/t_300
err_rate_4 = np.sqrt(rate_4)/t_300

#  Logarithmen draufschmeißen
I_1 = np.log(I_0/rate_1)
I_2 = np.log(I_0/rate_2)
I_4 = np.log(I_0/rate_4)

err_I_1 = 1/I_0 * err_I_0 - 1/rate_1 * err_rate_1
err_I_2 = 1/I_0 * err_I_0 - 1/rate_2 * err_rate_2
err_I_4 = 1/I_0 * err_I_0 - 1/rate_4 * err_rate_4

# Gewichtungsmatrizen
W_1 = np.diag(1/err_I_1**2)
W_2 = np.diag(1/err_I_2**2)
W_4 = np.diag(1/err_I_4**2)

# Über kleinsteQuadrate mu berechnen
mu_1, err_mu_1 = kleinsteQuadrate(I_1, W=W_1, A=A)
mu_2, err_mu_2 = kleinsteQuadrate(I_2, W=W_2, A=A)
mu_4, err_mu_4 = kleinsteQuadrate(I_4, W=W_4, A=A)

print(
    '''
~~~ Absorptionskoeffizienten der verschiedenen Würfel ~~~
Würfel 1:
-----------------------------------------------
Werte = {}
Fehler = {}
Würfel 2:
-----------------------------------------------
Werte = {}
Fehler = {}
Würfel 3:
-----------------------------------------------
Werte = {}
Fehler = {}
'''.format(mu_1, err_mu_1, mu_2, err_mu_2, mu_4, err_mu_4)
)

mu_1_mit = unumpy.uarray(mu_1, err_mu_1)
mu_2_mit = unumpy.uarray(mu_2, err_mu_2)
mu_4_mit = unumpy.uarray(mu_4, err_mu_4)

mu_1_mittel = np.mean(mu_1_mit)
mu_2_mittel = np.mean(mu_2_mit)
print("Würfel 1: ", mu_1_mittel)
print("Würfel 2: ", mu_2_mittel)
lol = np.linspace(1, 9, 9)

# Literaturwerte
rho = np.array([2.699, 11.35, 7.874, 8.284, 1.405])
sigma = np.array([7.466*10**(-2), 1.101*10**(-1), 7.346*10**(-2),
                  7.514*10**(-2), 0.086])

mu = rho * sigma
print(mu)

# Abweichungen
# Für jeden Elementarwürfel und jedes Element die relative Abweichung bestimmen
abweichungen = np.array([])
for i in range(len(mu_4_mit)):
    for j in range(len(mu)):
        x = mu_4_mit[i] / mu[j]
        if x > 1:
            x = x-1
            abweichungen = np.append(abweichungen, x)
        else:
            x = 1-x
            abweichungen = np.append(abweichungen, x)

# Nach Elementarwürfeln aufteilen und dabei das kleinste Element und dessen
# Position übergeben
abweichungen_1 = np.array([min(abweichungen[0:5])*100,
                           np.argmin(abweichungen[0:5])])
abweichungen_2 = np.array([min(abweichungen[5:10])*100,
                           np.argmin(abweichungen[5:10])])
abweichungen_3 = np.array([min(abweichungen[10:15])*100,
                           np.argmin(abweichungen[10:15])])
abweichungen_4 = np.array([min(abweichungen[15:20])*100,
                           np.argmin(abweichungen[15:20])])
abweichungen_5 = np.array([min(abweichungen[20:25])*100,
                           np.argmin(abweichungen[20:25])])
abweichungen_6 = np.array([min(abweichungen[25:30])*100,
                           np.argmin(abweichungen[25:30])])
abweichungen_7 = np.array([min(abweichungen[30:35])*100,
                           np.argmin(abweichungen[30:35])])
abweichungen_8 = np.array([min(abweichungen[35:40])*100,
                           np.argmin(abweichungen[35:40])])
abweichungen_9 = np.array([min(abweichungen[40:45])*100,
                           np.argmin(abweichungen[40:45])])

print("Abweichung erster Würfel: ", abweichungen_1)
print("Abweichung zweiter Würfel: ", abweichungen_2)
print("Abweichung dritter Würfel: ", abweichungen_3)
print("Abweichung vierter Würfel: ", abweichungen_4)
print("Abweichung fünfter Würfel: ", abweichungen_5)
print("Abweichung sechster Würfel: ", abweichungen_6)
print("Abweichung siebter Würfel: ", abweichungen_7)
print("Abweichung achter Würfel: ", abweichungen_8)
print("Abweichung neunter Würfel: ", abweichungen_9)

print((1 - 0.202/mu_2_mittel)*100)
# abweichungen_kleinste = np.array([abweichungen_1, abweichungen_2,
#                                  abweichungen_3, abweichungen_4,
#                                  abweichungen_5, abweichungen_6,
#                                  abweichungen_7, abweichungen_8,
#                                  abweichungen_9])

# Tabellen
np.savetxt('Tabellen/absorp_w1.txt',
           np.column_stack([lol, mu_1*10**3, err_mu_1*10**3]),
           delimiter=' \pm ', newline=r' \\'+'\n', fmt="%.2f")

np.savetxt('Tabellen/absorp_w2.txt',
           np.column_stack([lol, mu_2*10**3,
                            err_mu_2*10**3]),
           delimiter=' \pm ', newline=r' \\'+'\n', fmt="%.2f")

np.savetxt('Tabellen/absorp_w4.txt',
           np.column_stack([lol, mu_4*10**3, err_mu_4*10**3]),
           delimiter=' \pm ', newline=r' \\'+'\n', fmt="%.2f")

np.savetxt('Tabellen/mu.txt',
           np.column_stack([sigma, rho, mu]),
           delimiter=' & ', newline=r' \\'+'\n', fmt="%.3f")

print('Alles ausgeführt!')
