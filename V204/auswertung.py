import numpy as np
import uncertainties.unumpy as unp
from scipy import stats
from uncertainties import ufloat

#Auswertung für Messing
Brass_Anah_cm, Brass_Afern_cm, Brass_dt_cm = np.genfromtxt('data_brass', unpack=True)

Brass_ln = np.log((Brass_Anah_cm * 10 * 0.44)/(Brass_Afern_cm * 10 * 0.44))
Brass_ln_fehler = ufloat(np.mean(Brass_ln), stats.sem(Brass_ln))

Brass_dt_s = ufloat(np.mean(Brass_dt_cm * 10 * 4.35), stats.sem(Brass_dt_cm * 10 * 4.35))

dx = 0.03
rho_brass = 8520
c_brass = 385

print(Brass_ln_fehler)
print(Brass_dt_s)

print( ( rho_brass * c_brass * (dx)**2 ) / ( 2 * Brass_dt_s * Brass_ln_fehler ) )
#print(unp.log(Alu_Anah_K / Alu_Afern_K))
print("")

#np.savetxt('Brass.txt', np.column_stack([Brass_Anah_cm * 10 * 0.44, Brass_Afern_cm * 10 * 0.44, Brass_dt_cm * 10 * 4.35, Brass_ln]), fmt="%.2f")

#Auswertung für Aluminium
Alu_Anah_cm, Alu_Afern_cm, Alu_dt_cm = np.genfromtxt('data_alu', unpack=True)

Alu_ln = np.log((Alu_Anah_cm * 10 * 0.37)/(Alu_Afern_cm * 10 * 0.37))
Alu_ln_fehler = ufloat(np.mean(Alu_ln), stats.sem(Alu_ln))

Alu_dt_s = ufloat(np.mean(Alu_dt_cm * 10 * 4.44), stats.sem(Alu_dt_cm * 10 * 4.44))

dx = 0.03
rho_alu = 2800
c_alu = 830

print(Alu_ln_fehler)
print(Alu_dt_s)

print( ( rho_alu * c_alu * (dx)**2 ) / ( 2 * Alu_dt_s * Alu_ln_fehler ) )
#print(unp.log(Alu_Anah_K / Alu_Afern_K))
print("")

#np.savetxt('Alu.txt', np.column_stack([Alu_Anah_cm * 10 * 0.37, Alu_Afern_cm * 10 * 0.37, Alu_dt_cm * 10 * 4.44, Alu_ln]), fmt="%.2f")

#Auswertung für Edelstahl
Ssteel_Anah_cm, Ssteel_Afern_cm, Ssteel_dt_cm = np.genfromtxt('data_ssteel', unpack=True)

Ssteel_ln = np.log((Ssteel_Anah_cm * 10 * 0.11)/(Ssteel_Afern_cm * 10 * 0.11))
Ssteel_ln_fehler = ufloat(np.mean(Ssteel_ln), stats.sem(Ssteel_ln))

Ssteel_dt_s = ufloat(np.mean(Ssteel_dt_cm * 10 * 9.10), stats.sem(Ssteel_dt_cm * 10 * 9.10))

dx = 0.03
rho_ssteel = 8000
c_ssteel = 400

print(Ssteel_ln_fehler)
print(Ssteel_dt_s)

print( ( rho_ssteel * c_ssteel * (dx)**2 ) / ( 2 * Ssteel_dt_s * Ssteel_ln_fehler ) )


#np.savetxt('Ssteel.txt', np.column_stack([Ssteel_Anah_cm * 10 * 0.11, Ssteel_Afern_cm * 10 * 0.11, Ssteel_dt_cm * 10 * 9.10, Ssteel_ln ]), fmt="%.2f")
