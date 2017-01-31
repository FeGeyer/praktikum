import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
#from uncertainties import uarray
from scipy.optimize import curve_fit
from scipy.stats import stats
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 15

#Einlesen der Daten
Ck_a_b_in_nF, Bäuche, fehler_bäuche, v_plus_in_kHz, v_minus_in_kHz, fehler_frequenzen = np.genfromtxt('messdaten_a_b.txt', unpack=True)
Ck_c_in_nF, Delta_t_1_in_ms, Delta_t_2_in_ms, fehler_zeiten = np.genfromtxt('messdaten_c.txt', unpack=True)
L_in_mH = 23.954
C_in_nF = 0.7932
Csp_in_nF = 0.028
f_start_in_kHz = 20.06
f_end_in_kHz = 68.97
f_resonanz_gemessen_in_kHz = 35.65
Fehler_Ck = 0.003
f_res_gemessen = ufloat(35.65, 0.02)

#Messgrößen mit Fehlern
Bäuche = unp.uarray(Bäuche, fehler_bäuche)
v_plus_in_kHz = unp.uarray(v_plus_in_kHz, fehler_frequenzen)
v_minus_in_kHz = unp.uarray(v_minus_in_kHz, fehler_frequenzen)
Delta_t_1_in_ms = unp.uarray(Delta_t_1_in_ms, fehler_zeiten)
Delta_t_2_in_ms = unp.uarray(Delta_t_2_in_ms, fehler_zeiten)

#Kapazitäten mit Fehlern (gerundet auf 3 signifikante Stellen)
Ck_a_b_Fehler_nF = unp.uarray(Ck_a_b_in_nF, np.round(Ck_a_b_in_nF * Fehler_Ck, 3))
Ck_c_Fehler_nF = unp.uarray(Ck_c_in_nF, np.round(Ck_c_in_nF * Fehler_Ck, 3))

#Werte in SI
Ckab = Ck_a_b_Fehler_nF * 10**(-9)
Ckc = Ck_c_Fehler_nF * 10**(-9)
L = L_in_mH * 10**(-3)
C = C_in_nF * 10**(-9)
Csp = Csp_in_nF * 10**(-9)

#Bestimmen der Resonanzfrequenz
f_res = 1/(2 * np.pi * np.sqrt(L_in_mH * 10**(-3) * C_in_nF * 10**(-9) ))
#print(f_res)
#print(100*(1-(f_res_gemessen*10**(3) / f_res)))

#Bestimmen der rechnerischen Frequenzen nu+ und nu-
nu_plus = 1/(2 * np.pi * np.sqrt(L * (C+ Csp)) )
nu_minus = 1/(2 * np.pi * unp.sqrt( L * ( Csp + (1/C + 2/Ckab)**(-1) ) ) )

#Berechnen der rechnerischen Anzahl der Maxima aus nu+ und nu-
Bäuche_rechnerisch = (nu_plus + nu_minus)  / ((nu_minus - nu_plus) )
print(nu_plus)
print(nu_minus)
print(Bäuche_rechnerisch)
rel_Bäuche = (1 - Bäuche / Bäuche_rechnerisch) * 100


#Berechnen der relativen Abweichungen zwischen gemessener Fundamentalschwingung und rechnerischen
rel_nu_plus = 1 - (v_plus_in_kHz * 10**(3)) / nu_plus
rel_nu_minus = 1 - nu_minus / (v_minus_in_kHz * 10**(3))
print(rel_nu_plus * 100)
print(rel_nu_minus * 100)
np.savetxt('Tabelle_a_b.txt', np.column_stack([unp.nominal_values(Ck_a_b_in_nF), unp.nominal_values(v_minus_in_kHz), unp.nominal_values(nu_minus) * 10**(-3),
 unp.nominal_values(rel_nu_minus * 100)]), fmt="%.2f" )
np.savetxt('Tabelle_a_b_2.txt', np.column_stack([unp.nominal_values(Ck_a_b_in_nF), unp.nominal_values(Bäuche), unp.nominal_values(Bäuche_rechnerisch),
 unp.std_devs(Bäuche_rechnerisch), unp.nominal_values(rel_Bäuche)]), fmt="%.2f" )
#Brechnen der Fundamentalfreuqenzen aus der Sweep-Methode
rech_sweep_nu_plus = 1/(2 * np.pi * np.sqrt(L * (C+ Csp)) )
rech_sweep_nu_minus = 1/(2 * np.pi * unp.sqrt( L * ( Csp + (1/C + 2/Ckc)**(-1) ) ) )
sweep_nu_plus = f_start_in_kHz * 10**(3) + ( (f_end_in_kHz - f_start_in_kHz) * 10**(3) * Delta_t_1_in_ms * 10**(-3) )
sweep_nu_minus = f_start_in_kHz * 10**(3) + ( (f_end_in_kHz - f_start_in_kHz) * 10**(3) * Delta_t_2_in_ms * 10**(-3) )


rel_sweep_nu_plus = 1 - (sweep_nu_plus / rech_sweep_nu_plus)
rel_sweep_nu_minus = 1 - (rech_sweep_nu_minus / sweep_nu_minus)

print(sweep_nu_plus)
print(rel_sweep_nu_plus)
print(rel_sweep_nu_minus)

np.savetxt('Tabelle_c.txt', np.column_stack([unp.nominal_values(Ck_c_in_nF), unp.nominal_values(Delta_t_2_in_ms), unp.nominal_values(sweep_nu_minus)*10**(-3),
unp.nominal_values(rech_sweep_nu_minus) * 10**(-3),
unp.std_devs(rech_sweep_nu_minus) * 10 ** (-3), unp.nominal_values(rel_sweep_nu_minus * 100)]), fmt="%.2f" )

print(sweep_nu_minus*10**(-3))

#Plot
t = np.linspace(0.9 * 10 ** (-9), 12.1 * 10 ** (-9), 100)
plt.figure(1)
plt.xlabel(r"$C_k / \mathrm{nF}$")
plt.ylabel(r"$\nu_- / \mathrm{kHz}$")
plt.grid(True, which="both")
plt.ylim(36,57)
plt.xlim(0.9, 12.1)
plt.plot(t * 10 ** 9, unp.nominal_values(1/(2 * np.pi * unp.sqrt(L * (Csp + (1/C + 2/t) ** (-1) ) ) )) * 10 ** (-3), 'r-', label='Theoriekurve')
plt.plot(t, t**2, "r+")
plt.plot(unp.nominal_values(Ck_a_b_in_nF), unp.nominal_values(v_minus_in_kHz), 'g+', label='Lissajous-Methode')
plt.plot(unp.nominal_values(Ck_c_in_nF), unp.nominal_values(sweep_nu_minus) * 10**(-3), 'b+', label='Sweep-Methode')
plt.legend(loc="best")
plt.tight_layout
plt.savefig('Methoden.pdf')
