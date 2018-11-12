from uncertainties import ufloat
from uncertainties import unumpy
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

# Import GreyValue Data. Keys: Pixel, GreyValue
Metall = pd.read_csv(filepath_or_buffer="Bilder/Metall.csv")

# Convert from Pixel to centimetre, distance mesuered to be 16.5 cm
Metall.Distance = Metall.Pixel*(18/(len(Metall.Pixel)+6))
Metall.GreyValue = -Metall.GreyValue

Metall.GreyValue[Metall.GreyValue < -18] = -18

# Find peaks
MetallPeaks, _ = find_peaks(x=Metall.GreyValue, prominence=2.5)
GreyValuePeaks = Metall.GreyValue[MetallPeaks]

# Get Distance from peaks
MetallPeaks = MetallPeaks*(18/(len(Metall.Pixel)+6))

plt.figure()
plt.plot(Metall.Distance, Metall.GreyValue, ls='--', color='blue')
plt.plot(MetallPeaks, GreyValuePeaks, color='black', ls='', marker='o')
plt.xlabel(r"$r / \mathrm{cm}$")
plt.ylabel('inverser Grauwert')
plt.xlim(0, 18)
plt.show()

# Distance is equal to 180° -> 1cm equals 10°
# Bragg: 2dsin(theta)=n*lambda
# lambda = 1.54093A
# Angles have to be reverted, cause they have to be mesured acording to
# the MP-vector

lam = 1.54093*10**(-10)
R = 57.3*10**(-3)

PeakAngle = MetallPeaks/(2*R)
PeakAngle = np.abs(PeakAngle-180)

print(PeakAngle)

# The first reflexes are due to the existence of to K-alpha lines, so correct for it

# theta = np.arctan((PeakAngle[0]-PeakAngle[1])*(1.5417*10**(-10))/(1.54478*10**(-10)-lam))
theta = (PeakAngle[0]-PeakAngle[1])/2
PeakAngle = np.insert(PeakAngle, 0, PeakAngle[0]-theta)
PeakAngle = np.delete(PeakAngle, [1, 2])

print(PeakAngle)
print(theta)

d = lam/(2*np.sin(PeakAngle/2))

# Possible Millerindizes for fcc and bcc
############
# 1: bcc: h+k+l gerade
# 2: fcc: h,k,l alle gerade oder ungerade
############

# nenner1 = np.sqrt(h1**2 + k1**2 + l1**2)
# nenner2 = np.sqrt(h2**2 + k2**2 + l2**2)
#
# a1 = d*nenner1
# a2 = d*nenner2

print(a1)
print(a2)
