import numpy as np
import uncertainties.unumpy as unp
from uncertainties import ufloat
import scipy.constants as const
from scipy.constants import physical_constants

m = ufloat(4.6913, 0)
I = ufloat(8.3333e-10, 0)
a = ufloat(0.0101695151722, 4.74488349932e-05)
g = const.g

E = (m * g)/(48 * I * a)

print(E)
print(const.g)
