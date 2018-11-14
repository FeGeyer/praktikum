from scipy.signal import find_peaks
from scipy.stats import sem
import latticeconstants as lat
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

# Import GreyValue Data. Keys: Pixel, GreyValue
Metall = pd.read_csv(filepath_or_buffer="Bilder/Metall.csv")

# Convert from Pixel to centimetre, distance mesuered to be 16.5 cm
Metall.Distance = Metall.Pixel * (18 / (len(Metall.Pixel) + 6))
Metall.GreyValue = -Metall.GreyValue

Metall.GreyValue[Metall.GreyValue < -18] = -18

# Find peaks
MetallPeaks, _ = find_peaks(x=Metall.GreyValue, prominence=2.5)
GreyValuePeaks = Metall.GreyValue[MetallPeaks]

# Get Distance from peaks
MetallPeaks = MetallPeaks * (18 / (len(Metall.Pixel) + 6))

# plt.figure()
# plt.plot(Metall.Distance, Metall.GreyValue, ls='--', color='blue')
# plt.plot(MetallPeaks, GreyValuePeaks, color='black', ls='', marker='o')
# plt.xlabel(r"$r / \mathrm{cm}$")
# plt.ylabel('inverser Grauwert')
# plt.xlim(0, 18)
# plt.show()

# Distance is equal to 180° -> 1cm equals 10°
# Bragg: 2dsin(theta)=n*lambda
# lambda = 1.54093A
# Angles have to be reverted, cause they have to be mesured acording to
# the MP-vector

lam = 1.54093 * 10**(-10)
R = 57.3 * 10**(-3)

PeakAngle = MetallPeaks * 10
PeakAngle = np.abs(180 - PeakAngle)
PeakAngle = np.sort(PeakAngle)

# The first two reflexes are due to the existence of two K-alpha lines,
# so correct for it
theta = (PeakAngle[0] - PeakAngle[1]) / 2
PeakAngle = np.insert(PeakAngle, 0, PeakAngle[0] - theta)
PeakAngle = np.delete(PeakAngle, [1, 2])

print(PeakAngle)
print(theta)

# convert the angles to interplanar distance according to braggs law
d = lam / (2 * np.sin(0.5 * PeakAngle * np.pi / 180))

# Compute possible millerindizes for fcc and bcc
h_bcc, k_bcc, l_bcc = lat.bcc(7)
h_fcc, k_fcc, l_fcc = lat.fcc(7)

n_fcc = np.sqrt(h_bcc**2 + k_bcc**2 + l_bcc**2)
n_bcc = np.sqrt(h_fcc**2 + k_fcc**2 + l_fcc**2)


def find_best_fit(evaluate_value, test_array, i):

    '''
    evaluate_value: float, value, for witch the best fitting value for n
    has to be evaluated.
    test_array: array of floats, possible values for evaluate_value
    i: int, index of currently evaluated PeakAngle

    Function to find the best fitting value, therfor the value, where the
    euclidian distance between evaluate_value and the evaluated test_array
    value is minimal.

    Returns best fit and distance to evaluated value (residuum res)
    '''

    res = d[i] / np.abs(test_array - evaluate_value)
    best_fit = test_array[d[i] / np.abs(test_array - evaluate_value) ==
                          d[i] / np.abs(test_array - evaluate_value).min()]
    best_res = res[d[i] / np.abs(test_array - evaluate_value) ==
                   d[i] / np.abs(test_array - evaluate_value).min()]
    return best_fit[0], best_res[0]

# Combining Bragg and the interplanar distance latticeconstant relation gives
# sin(theta_2)/sin(theta_1) = n_1/n_2 = c with n = np.sqrt(h**2+k**2+l**2)

# initialize array for c
c_bcc = np.zeros([len(PeakAngle) - 1])

# compute c values
for angle_index in range(len(PeakAngle) - 1):
    c_bcc[angle_index] = (np.sin(PeakAngle[angle_index + 1] * np.pi / 180) /
                          np.sin(PeakAngle[angle_index] * np.pi / 180))

# initialize field, where the n_values for each initial guess of n are
# computeted iterative and stored
n_search_field = np.zeros([len(PeakAngle), len(n_bcc)])

# set initial guess on all possible n
n_search_field[0, :] = n_bcc

# iterativly compute following n
for row in range(len(PeakAngle) - 1):
    for column in range(len(n_bcc)):
        n_search_field[row + 1, column] = (n_search_field[row, column] /
                                           c_bcc[row])

# seafty first, even if negative values should not happen
# n_search_field = np.abs(n_search_field)
print(c_bcc)

# initialize fields for best fits and residdums, setting value of the initial
# guess according to the formally used values
n_best_field = np.zeros([len(PeakAngle), len(n_bcc)])
n_best_field[0, :] = n_bcc
n_res_field = np.zeros([len(PeakAngle), len(n_bcc)])
n_res_field[0, :] = 0

# compute best fits and risiduums, deltining formally used ns
for column in range(len(n_bcc)):
    n_bcc_help = n_bcc
    help_list = n_bcc_help.tolist()
    index = help_list.index(n_best_field[0, column])
    n_bcc_help = np.delete(n_bcc_help, [index])
    for row in range(len(PeakAngle) - 1):
        (n_best_field[row + 1, column],
            n_res_field[row + 1, column]) = find_best_fit(n_search_field[
                                                          row + 1, column],
                                                          n_bcc_help,
                                                          row + 1)
        help_list = n_bcc_help.tolist()
        index = help_list.index(n_best_field[row + 1, column])
        n_bcc_help = np.delete(n_bcc_help, [index])
# compute mean of residuums per column and find minimal value
res_sum = np.mean(n_res_field, axis=0)
res_sem = sem(n_res_field, axis=0)
res_sum_list = res_sum.tolist()
best_index = res_sum_list.index(res_sum.min())
best_sum = res_sum[best_index]
best_sem = res_sem[best_index]

print(n_best_field[:, best_index]**2)
print(d / n_best_field[:, best_index])
print(best_sum, best_sem)
