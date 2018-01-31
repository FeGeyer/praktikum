import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit

D_theta, R1, R2, R3, R4, R5 = np.genfromtxt('kontrast.txt', unpack=True)

T = 1 * 10**-3
lam = 633 *10**-9
tetha_0 =
Winkel = Winkela * (np.pi)/180
Winkel2 = Winkelb * (np.pi)/180
Winkel3 = Winkelc * (np.pi)/180


n1 = (- T * 2* tetha_0 * Winkel) / (Anzahl * lam  - T*(2*tetha_0 *Winkel))
