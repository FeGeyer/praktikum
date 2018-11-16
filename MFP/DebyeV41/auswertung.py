from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.stats import sem
from uncertainties import ufloat
import matplotlib.pyplot as plt
'''
Import millerindizes as miller from Auswertung/millerindizes.py with
importlib.util. Important to note is, that in the third line of this piece of
code the name of the module is definded for further use in the code.
'''
import importlib.util
spec = importlib.util.spec_from_file_location("millerindizes", "Auswertung/millerindizes.py")
miller = importlib.util.module_from_spec(spec)
spec.loader.exec_module(miller)

import numpy as np
import pandas as pd
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

# Turn of pandas warnings... be carefull
pd.options.mode.chained_assignment = None


def linear(x, m, n):

    '''
    x: array of floats, values, for wich the duction has to be evaluated.
    m, n: floats, Parameter for the function.

    Simple linear function.

    Returns fuction values y=f(x).
    '''

    return m*x+n


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
    best_fit = test_array[d[i] * np.abs(test_array - evaluate_value) ==
                          d[i] * np.abs(test_array - evaluate_value).min()]
    best_res = res[d[i] * np.abs(test_array - evaluate_value) ==
                   d[i] * np.abs(test_array - evaluate_value).min()]
    return best_fit[0], best_res[0]


def find_lattice_constants(d, lattice, max_value):

    '''
    d: array of floats, array of the interplanar distances
    lattice
    lattice: string, assumed lattice-type. Supported: bcc and fcc
    max_value: int, maximum value for h, k and l respectivly

    Funktion to find the best fitting millier indizes (h, k, l) and compute the
    lattice-constant a.

    Returns best fitting n = sqrt(h**2 + k**2 + l**2), the a for any n and
    mean a and its error.
    '''

    # Compute possible millerindizes for given lattice and max_value
    if lattice == "bcc":
        h, k, l = miller.bcc(max_value)
    elif lattice == "fcc":
        h, k, l = miller.fcc(max_value)
    else:
        print("No supported lattice-type given")
        return

    # Compute denominator of latticeconstant formular
    n = np.sqrt(h**2 + k**2 + l**2)

    # Combining Bragg and the interplanar distance
    # latticeconstant relation gives
    # sin(theta_2)/sin(theta_1) = n_1/n_2 = c with n = np.sqrt(h**2+k**2+l**2)

    # initialize array for c
    c_bcc = np.zeros([len(PeakAngle) - 1])

    # compute c values
    for angle_index in range(len(PeakAngle) - 1):
        c_bcc[angle_index] = (np.sin(PeakAngle[angle_index + 1] * np.pi /
                              180) / np.sin(PeakAngle[angle_index] * np.pi /
                              180))

    # initialize field, where the n_values for each initial guess of n are
    # computeted iterative and stored
    n_search_field = np.zeros([len(PeakAngle), len(n)])

    # set initial guess on all possible n
    n_search_field[0, :] = n

    # iterativly compute following n
    for row in range(len(PeakAngle) - 1):
        for column in range(len(n)):
            n_search_field[row + 1, column] = (n_search_field[row, column] /
                                               c_bcc[row])

    # seafty first, even if negative values should not happen
    # n_search_field = np.abs(n_search_field)
    print(c_bcc)

    # initialize fields for best fits and residdums,
    # setting value of the initial
    # guess according to the formally used values
    n_best_field = np.zeros([len(PeakAngle), len(n)])
    n_best_field[0, :] = n
    n_res_field = np.zeros([len(PeakAngle), len(n)])
    n_res_field[0, :] = 0

    # compute best fits and risiduums, deltining formally used ns
    for column in range(len(n)):
        n_help = n
        help_list = n_help.tolist()
        index = help_list.index(n_best_field[0, column])
        n_help = np.delete(n_help, [index])
        for row in range(len(PeakAngle) - 1):
            (n_best_field[row + 1, column],
             n_res_field[row + 1, column]) = find_best_fit(n_search_field[
                                                           row + 1, column],
                                                           n_help,
                                                           row + 1)
            help_list = n_help.tolist()
            index = help_list.index(n_best_field[row + 1, column])
            n_help = np.delete(n_help, [index])
    # compute mean of residuums per column and find minimal value
    res_sum = np.mean(n_res_field, axis=0)
    res_sem = sem(n_res_field, axis=0)
    res_sum_list = res_sum.tolist()
    best_index = res_sum_list.index(res_sum.min())
    best_sum = res_sum[best_index]
    best_sem = res_sem[best_index]

    n_bestimmt = n_best_field[:, best_index]
    a_bestimmt = d * n_best_field[:, best_index]

    return n_bestimmt, a_bestimmt, best_sum, best_sem


if __name__ == '__main__':
    # Print a welcome text.
    print("Welcome pesant or creator, whatever holds true. If you are the "
          "author of these script, you can scip this lines. "
          "If not, you may find some muy bien importante informationes here: "
          "First of all, you will need to produce GreyValue "
          "Data for this script. "
          "We recommend the use of fiji (ImageJ), where the "
          "Analysis -> Plot Profile "
          "function gives yuo a GreyValue Plot whose data can "
          "be saved as .csv file. "
          "The best results are produced, if the films "
          "are scanned. It might be "
          "necessary to leave the scanner open and expose the film to light "
          "(e.g. like a flashlight). Second to that, you "
          "will need to adjust some "
          "lines in this script. They are marked with a blockcomment.")
    print("------------------------------------------------------------------")

    '''
    main-routine of the analysis. In here, data is beeing read in, peaks are
    computed and the latticeconstants are computed.
    For this, the functions above are called.
    '''
    # Metall-Probe
    # Import GreyValue Data. Keys: Pixel, GreyValue
    Metall = pd.read_csv(filepath_or_buffer="Auswertung/Bilder/Metall.csv")
    Salz = pd.read_csv(filepath_or_buffer="Auswertung/Bilder/Salz.csv")

    # Subtract underground
    Salz.GreyValue = -Salz.GreyValue
    Metall.GreyValue = -Metall.GreyValue

    '''
    Attention:
    In the following rows, it is indispensable to adjust the cut-of value
    inside the boolean statement of the mask to your needs.
    '''

    Salz.GreyValue[Salz.GreyValue < -175] = -175
    Metall.GreyValue[Metall.GreyValue < -18] = -18

    for idx, Probe in enumerate([Metall, Salz]):

        '''
        Attention:
        In the following rows, the summands 6 (this happens twice!)
        is used to correct for a black
        marker, placed inside the center of the punchholes of the film to
        simplify the selection of a alaysis-window in ImageJ .
        Since they were colored black, the could not be taken in account
        for the grayvalue messurement. To evoid an error through this, it
        is corrected by adding the thereby lost 6 pixels in the conversion.
        '''
        # Convert from Pixel to centimetre, distance mesuered to be 18 cm
        print("Just ignore this warning, it seems to be useless:")
        Probe.Distance = Probe.Pixel * (18 / (len(Probe.Pixel) + 6))
        print("--------------------------------------------------------------")

        # Find peaks
        ProbePeaks, _ = find_peaks(x=Probe.GreyValue, prominence=2.5)
        GreyValuePeaks = Probe.GreyValue[ProbePeaks]

        # Get Distance from peaks
        ProbePeaks = ProbePeaks * (18 / (len(Probe.Pixel) + 6))

        # Plot peaks
        plt.figure()
        plt.plot(Probe.Distance, Probe.GreyValue, ls='--', color='blue')
        plt.plot(ProbePeaks, GreyValuePeaks, color='black',
                 ls='', marker='o')
        plt.xlabel(r"$r / \mathrm{cm}$")
        plt.ylabel('inverser Grauwert')
        plt.xlim(0, 18)
        plt.show()

        # Distance is equal to 180° -> 1cm equals 10°
        # Bragg: 2dsin(theta)=n*lambda
        # lambda = 1.54093A
        # Angles have to be reverted, cause they have to be mesured acording to
        # the MP-vector

        lam = 1.54093 * 10**(-10)
        R = 57.3 * 10**(-3)

        PeakAngle = ProbePeaks * 10
        '''
        Attention:
        We had to invert the film in the next rows.
        Check if you need to do this.
        '''
        PeakAngle = np.abs(180 - PeakAngle)
        PeakAngle = np.sort(PeakAngle)

        # convert the angles to interplanar distance according to braggs law
        d = lam / (2 * np.sin(0.5 * PeakAngle * np.pi / 180))

        print("bcc c=sqrt(h**2+k**2+l**2):")
        bcc_n, bcc_a, bcc_mean, bcc_sem = find_lattice_constants(d, 'bcc', 7)
        print("fcc c=sqrt(h**2+k**2+l**2):")
        fcc_n, fcc_a, fcc_mean, fcc_sem = find_lattice_constants(d, 'fcc', 7)

        # Compute best_fit for a thrugh linear regression
        bcc_params, cov = curve_fit(linear, np.cos(0.5 * PeakAngle * np.pi /
                                    180)**2, bcc_a)
        bcc_errors = np.sqrt(np.diag(cov))
        m_bcc = ufloat(bcc_params[0], bcc_errors[0])
        n_bcc = ufloat(bcc_params[1], bcc_errors[1])

        fcc_params, cov = curve_fit(linear, np.cos(0.5 * PeakAngle * np.pi /
                                    180)**2, fcc_a)
        fcc_errors = np.sqrt(np.diag(cov))
        m_fcc = ufloat(fcc_params[0], fcc_errors[0])
        n_fcc = ufloat(fcc_params[1], fcc_errors[1])

        print("bcc best fit:")
        print("m = ", m_bcc, ", n = ", n_bcc)
        print("fcc best fit:")
        print("m = ", m_fcc, ", n = ", n_fcc)

        x_range = np.linspace(0, 90, 1000)
        x_range = np.cos(x_range * np.pi / 180)**2

        # Plot peaks
        plt.figure()
        plt.plot(np.cos(PeakAngle * 0.5 * np.pi / 180)**2, bcc_a * 10**(12),
                 marker='x', color='red', ls='')
        plt.plot(x_range, linear(x_range, *bcc_params) * 10**(12),
                 ls='-', color='red', label='Hypothese: bcc-Gitter')
        plt.plot(np.cos(PeakAngle * 0.5 * np.pi / 180)**2, fcc_a * 10**(12),
                 marker='x', color='blue', ls='')
        plt.plot(x_range, linear(x_range, *fcc_params) * 10**(12),
                 ls='-', color='blue', label='Hypothese: fcc-Gitter')
        plt.xlabel(r"$\cos{(\phi)}^{2}$")
        plt.ylabel(r'Berechnete Gitterkonstante$ / \mathrm{pm}$')
        plt.legend(loc="best")
        plt.xlim(0, 1)
        plt.tight_layout
        plt.show()

print("------------------------------------------------------------------")
print('Thats all folks!')
