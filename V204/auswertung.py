import numpy as np
import uncertainties.unumpy as unp
from scipy import stats
from uncertainties import ufloat

#Auswertung für Messing
Brass_Anah_cm, Brass_Afern_cm, Brass_dt_cm = np.genfromtxt('data_brass', unpack=True)

Brass_Anah_K = ufloat(np.mean(Brass_Anah_cm * 10 * 0.44), stats.sem(Brass_Anah_cm * 10 * 0.44))
Brass_Afern_K = ufloat(np.mean(Brass_Afern_cm * 10 * 0.44), stats.sem(Brass_Afern_cm * 10 * 0.44))
Brass_dt_s = ufloat(np.mean(Brass_dt_cm * 10 * 4.35), stats.sem(Brass_dt_cm * 10 * 4.35))

dx = 0.03
rho_brass = 8520
c_brass = 385

print(Brass_Anah_K)
print(Brass_Afern_K)
print(Brass_dt_s)

print( ( rho_brass * c_brass * (dx)**2 ) / ( 2 * Brass_dt_s * unp.log(Brass_Anah_K / Brass_Afern_K) ) )
#print(unp.log(Alu_Anah_K / Alu_Afern_K))
print("")

#Auswertung für Aluminium
Alu_Anah_cm, Alu_Afern_cm, Alu_dt_cm = np.genfromtxt('data_alu', unpack=True)

Alu_Anah_K = ufloat(np.mean(Alu_Anah_cm * 10 * 0.37), stats.sem(Alu_Anah_cm * 10 * 0.37))
Alu_Afern_K = ufloat(np.mean(Alu_Afern_cm * 10 * 0.37), stats.sem(Alu_Afern_cm * 10 * 0.37))
Alu_dt_s = ufloat(np.mean(Alu_dt_cm * 10 * 4.44), stats.sem(Alu_dt_cm * 10 * 4.44))

dx = 0.03
rho_alu = 2800
c_alu = 830

print(Alu_Anah_K)
print(Alu_Afern_K)
print(Alu_dt_s)

print( ( rho_alu * c_alu * (dx)**2 ) / ( 2 * Alu_dt_s * unp.log(Alu_Anah_K / Alu_Afern_K) ) )
#print(unp.log(Alu_Anah_K / Alu_Afern_K))
print("")

#Auswertung für Edelstahl
Ssteel_Anah_cm, Ssteel_Afern_cm, Ssteel_dt_cm = np.genfromtxt('data_ssteel', unpack=True)

Ssteel_Anah_K = ufloat(np.mean(Ssteel_Anah_cm * 10 * 0.44), stats.sem(Ssteel_Anah_cm * 10 * 0.44))
Ssteel_Afern_K = ufloat(np.mean(Ssteel_Afern_cm * 10 * 0.44), stats.sem(Ssteel_Afern_cm * 10 * 0.44))
Ssteel_dt_s = ufloat(np.mean(Ssteel_dt_cm * 10 * 4.35), stats.sem(Ssteel_dt_cm * 10 * 4.35))

dx = 0.03
rho_ssteel = 8000
c_ssteel = 400

print(Ssteel_Anah_K)
print(Ssteel_Afern_K)
print(Ssteel_dt_s)

print( ( rho_ssteel * c_ssteel * (dx)**2 ) / ( 2 * Ssteel_dt_s * unp.log(Ssteel_Anah_K / Ssteel_Afern_K) ) )
