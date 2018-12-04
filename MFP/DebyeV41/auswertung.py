from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.stats import sem
from uncertainties import ufloat
import matplotlib.pyplot as plt
from matplotlib import pylab

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

    return m * x + n


def find_best_fit(evaluate_value, test_array, i):

    '''
    evaluate_value: float, value, for witch the best fitting value for n
    has to be evaluated.
    test_array: array of floats, possible values for evaluate_value
    i: int, index of currently evaluated PeakAngle

    Function to find the best fitting value, therfor the value, where the
    euclidian distance between evaluate_value and the evaluated test_array
    value is minimal. In other words, it just gives the lowest possible
    n, because the fundtion does not work like initially planned (due
    to a wrong implementation, to be honest). The reason it is not alterd
    is simply because it works.

    Returns best fit and distance to evaluated value (residuum res)
    '''

    # Initially, there should be a multiplication insted of a division.
    # But shit hits the fan if it is alterd that way.
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
    elif lattice == "Dia":
        h, k, l = miller.Dia(max_value)
    elif lattice == "CsCl":
        h, k, l = miller.CsCl(max_value)
    elif lattice == "ZnS":
        h, k, l = miller.ZnS(max_value)
    elif lattice == "F":
        h, k, l = miller.F(max_value)
    elif lattice == "NaCl":
        h, k, l = miller.NaCl(max_value)
    else:
        print("No supported lattice-type given")
        return

    # Compute denominator of latticeconstant formular
    n = np.sqrt(h**2 + k**2 + l**2)

    n = n[:len(d)]
    h = h[:len(d)]
    k = k[:len(d)]
    l = l[:len(d)]

    a = d * n
    return n, a


def find_hkl(lattice, n, max_value):

    '''
    lattice: string, assumed lattice-type. Supported: bcc and fcc
    n: float, sqrt(h**2 + k**2 + l**2)
    max_value: int, maximum value for h, k and l respectivly

    Funktion to find the corresponding millier indizes (h, k, l) to a n.

    Returns millerindizes.
    '''

    # Compute possible millerindizes for given lattice and max_value
    if lattice == "bcc":
        h, k, l = miller.bcc(max_value)
    elif lattice == "fcc":
        h, k, l = miller.fcc(max_value)
    elif lattice == "Dia":
        h, k, l = miller.Dia(max_value)
    elif lattice == "CsCl":
        h, k, l = miller.CsCl(max_value)
    elif lattice == "ZnS":
        h, k, l = miller.ZnS(max_value)
    elif lattice == "F":
        h, k, l = miller.F(max_value)
    elif lattice == "NaCl":
        h, k, l = miller.NaCl(max_value)
    else:
        print("No supported lattice-type given")
        return

    # Compute denominator of latticeconstant formular
    n_lattice = np.sqrt(h**2 + k**2 + l**2)
    indizes = np.empty([1])

    # Compare given n to possible n
    for j in range(len(n)):
        for i in range(len(n_lattice)):
            if n[j] == n_lattice[i]:
                # Set found n to 0 to avoid multiple identification of
                # the same h, k, l-tuple
                n_lattice[i] = 0
                indizes = np.append(indizes, i)
                break

    indizes = np.delete(indizes, [0])
    print(indizes)
    indizes = indizes.astype('int')

    # Identify indizes
    mask = np.zeros(len(h), dtype=bool)
    mask[indizes] = True

    h = h[mask]
    k = k[mask]
    l = l[mask]

    return h, k, l


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
    Metall.name = "Metall"
    Salt = pd.read_csv(filepath_or_buffer="Auswertung/Bilder/Salz.csv")
    Salt.name = "Salt"

    # Subtract underground
    Salt.GreyValue = -Salt.GreyValue
    Metall.GreyValue = -Metall.GreyValue

    '''
    Attention:
    In the following rows, it is indispensable to adjust the cut-of value
    inside the boolean statement of the mask to your needs.
    '''

    Salt.GreyValue[Salt.GreyValue < -175] = -175
    Metall.GreyValue[Metall.GreyValue < -18] = -18

    for idx, Probe in enumerate([Metall]):

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

        print(Probe.name, "probe")

        # Find peaks
        ProbePeaks, props = find_peaks(x=Probe.GreyValue, prominence=2.5)

        # Use only dark peaks
        LightPeaks = ProbePeaks[props['prominences'] <= 10]
        # ProbePeaks = ProbePeaks[props['prominences'] > 10]

        GreyValuePeaks = np.array(Probe.GreyValue[ProbePeaks])
        GreyValueLightPeaks = np.array(Probe.GreyValue[LightPeaks])

        # correct for dual peaks due to k_alpha_1 and k_alpha_2
        Corr_peak = (ProbePeaks[0] + ProbePeaks[1]) / 2
        Corr_Grey_Value = (GreyValuePeaks[0] + GreyValuePeaks[1]) / 2

        ProbePeaks = np.delete(ProbePeaks, [0, 1])
        GreyValuePeaks = np.delete(GreyValuePeaks, [0, 1])

        ProbePeaks = np.insert(ProbePeaks, [0], Corr_peak)
        GreyValuePeaks = np.insert(GreyValuePeaks, [0], Corr_Grey_Value)

        # Get Distance from peaks
        ProbePeaks = ProbePeaks * (18 / (len(Probe.Pixel) + 6))
        LightPeaks = LightPeaks * (18 / (len(Probe.Pixel) + 6))

        # Plot peaks
        plt.figure()
        plt.plot(Probe.Distance, Probe.GreyValue, ls='--', color='blue',
                 label="Grauwert")
        plt.plot(ProbePeaks, GreyValuePeaks, color='black',
                 ls='', marker='o', label="Erkannte Peaks")
        # plt.plot(LightPeaks, GreyValueLightPeaks, color='grey',
        #          ls='', marker='o', label="Erkannte, schwache Peaks")
        plt.xlabel(r"$r / \mathrm{cm}$")
        plt.ylabel('inverser Grauwert')
        plt.xlim(0, 18)
        plt.legend(loc="lower left")
        plt.tight_layout
        plt.savefig("Auswertung/Grafiken/" + Probe.name + "_Peaks.pdf")

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

        print("bcc n=sqrt(h**2+k**2+l**2):")
        bcc_n, bcc_a = find_lattice_constants(d, 'bcc', 7)
        print(bcc_n)
        print('bcc h, k, l:')
        bcc_h, bcc_k, bcc_l = find_hkl('bcc', bcc_n, 7)
        print(bcc_h, bcc_k, bcc_l)
        print("bcc a:")
        print(bcc_a)
        print("fcc n=sqrt(h**2+k**2+l**2):")
        fcc_n, fcc_a = find_lattice_constants(d, 'fcc', 7)
        print(fcc_n)
        print('fcc h, k, l:')
        fcc_h, fcc_k, fcc_l = find_hkl('fcc', fcc_n, 7)
        print(fcc_h, fcc_k, fcc_l)
        print("fcc a:")
        print(fcc_a)
        print("Dia n=sqrt(h**2+k**2+l**2):")
        Dia_n, Dia_a = find_lattice_constants(d, 'Dia', 7)
        print(Dia_n)
        print('Dia h, k, l:')
        Dia_h, Dia_k, Dia_l = find_hkl('Dia', Dia_n, 7)
        print(Dia_h, Dia_k, Dia_l)
        print("Dia a:")
        print(Dia_a)

        np.savetxt("Auswertung/Grafiken/" + Probe.name +
                   "_bcc_Tabelle.tex",
                   np.column_stack([
                                   PeakAngle,
                                   bcc_n**2,
                                   bcc_h,
                                   bcc_k,
                                   bcc_l,
                                   bcc_a * 10**(12),
                                   ]), delimiter=' & ', newline=r' \\' + '\n',
                   fmt='%.2f & %.0f & %.0f & %.0f & %.0f & %.2f')

        np.savetxt("Auswertung/Grafiken/" + Probe.name +
                   "_fcc_Tabelle.tex",
                   np.column_stack([
                                   PeakAngle,
                                   fcc_n**2,
                                   fcc_h,
                                   fcc_k,
                                   fcc_l,
                                   fcc_a * 10**(12),
                                   ]), delimiter=' & ', newline=r' \\' + '\n',
                   fmt='%.2f & %.0f & %.0f & %.0f & %.0f & %.2f')

        np.savetxt("Auswertung/Grafiken/" + Probe.name +
                   "_Dia_Tabelle.tex",
                   np.column_stack([
                                   PeakAngle,
                                   Dia_n**2,
                                   Dia_h,
                                   Dia_k,
                                   Dia_l,
                                   Dia_a * 10**(12),
                                   ]), delimiter=' & ', newline=r' \\' + '\n',
                   fmt='%.2f & %.0f & %.0f & %.0f & %.0f & %.2f')

        np.savetxt("Auswertung/Grafiken/" + Probe.name +
                   "_Tabelle.tex",
                   np.column_stack([
                                   PeakAngle,
                                   bcc_n**2,
                                   bcc_h,
                                   bcc_k,
                                   bcc_l,
                                   bcc_a * 10**(12),
                                   fcc_n**2,
                                   fcc_h,
                                   fcc_k,
                                   fcc_l,
                                   fcc_a * 10**(12),
                                   ]), delimiter=' & ', newline=r' \\' + '\n',
                   fmt='%.2f & %.0f & %.0f & %.0f & %.0f & %.2f & %.0f & %.0f & %.0f & %.0f & %.2f')

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

        Dia_params, cov = curve_fit(linear, np.cos(0.5 * PeakAngle * np.pi /
                                    180)**2, Dia_a)
        Dia_errors = np.sqrt(np.diag(cov))
        m_Dia = ufloat(Dia_params[0], Dia_errors[0])
        n_Dia = ufloat(Dia_params[1], Dia_errors[1])

        print("bcc best fit:")
        print("m = ", m_bcc, ", n = ", n_bcc)
        print("fcc best fit:")
        print("m = ", m_fcc, ", n = ", n_fcc)
        print("Dia best fit:")
        print("m = ", m_Dia, ", n = ", n_Dia)

        x_range = np.linspace(0, 90, 1000)
        x_range = np.cos(x_range * np.pi / 180)**2

        # Plot peaks
        plt.figure()
        plt.plot(np.cos(PeakAngle * 0.5 * np.pi / 180)**2, fcc_a * 10**(12),
                 marker='x', color='blue', ls='')
        plt.plot(x_range, linear(x_range, *fcc_params) * 10**(12),
                 ls='-', color='blue', label='Hypothese: fcc-Gitter')
        plt.xlabel(r"$\cos{(\phi)}^{2}$")
        plt.ylabel(r'Berechnete Gitterkonstante$ / \mathrm{pm}$')
        plt.legend(loc="best")
        plt.xlim(0, 1)
        plt.tight_layout
        plt.savefig("Auswertung/Grafiken/" +
                    Probe.name + "_fcc_Ausgleichsrechnung.pdf")

        plt.figure()
        plt.plot(np.cos(PeakAngle * 0.5 * np.pi / 180)**2, bcc_a * 10**(12),
                 marker='x', color='red', ls='')
        plt.plot(x_range, linear(x_range, *bcc_params) * 10**(12),
                 ls='-', color='red', label='Hypothese: bcc-Gitter')
        plt.xlabel(r"$\cos{(\phi)}^{2}$")
        plt.ylabel(r'Berechnete Gitterkonstante$ / \mathrm{pm}$')
        plt.legend(loc="best")
        plt.xlim(0, 1)
        plt.tight_layout
        plt.savefig("Auswertung/Grafiken/" +
                    Probe.name + "_bcc_Ausgleichsrechnung.pdf")
        plt.figure()
        plt.plot(np.cos(PeakAngle * 0.5 * np.pi / 180)**2, bcc_a * 10**(12),
                 marker='x', color='red', ls='')
        plt.plot(x_range, linear(x_range, *bcc_params) * 10**(12),
                 ls='-', color='red', label='Hypothese: bcc-Gitter')
        plt.plot(np.cos(PeakAngle * 0.5 * np.pi / 180)**2, fcc_a * 10**(12),
                 marker='x', color='blue', ls='')
        plt.plot(x_range, linear(x_range, *fcc_params) * 10**(12),
                 ls='-', color='blue', label='Hypothese: fcc-Gitter')
        # plt.plot(np.cos(PeakAngle * 0.5 * np.pi / 180)**2, Dia_a * 10**(12),
        #          marker='x', color='green', ls='')
        # plt.plot(x_range, linear(x_range, *Dia_params) * 10**(12),
        #          ls='-', color='green', label='Hypothese: Dia-Gitter')
        plt.xlabel(r"$\cos{(\phi)}^{2}$")
        plt.ylabel(r'Berechnete Gitterkonstante$ / \mathrm{pm}$')
        plt.legend(loc="best")
        plt.xlim(0, 1)
        plt.tight_layout
        plt.savefig("Auswertung/Grafiken/" +
                    Probe.name + "_Ausgleichsrechnung.pdf")

for idx, Probe in enumerate([Salt]):

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

    print(Probe.name, "probe")

    # Find peaks
    ProbePeaks, props = find_peaks(x=Probe.GreyValue, prominence=2.5)

    # Use only dark peaks
    LightPeaks = ProbePeaks[props['prominences'] <= 10]
    # ProbePeaks = ProbePeaks[props['prominences'] > 10]

    GreyValuePeaks = np.array(Probe.GreyValue[ProbePeaks])
    GreyValueLightPeaks = np.array(Probe.GreyValue[LightPeaks])

    # correct for dual peaks due to k_alpha_1 and k_alpha_2
    Corr_peak = (ProbePeaks[0] + ProbePeaks[1]) / 2
    Corr_Grey_Value = (GreyValuePeaks[0] + GreyValuePeaks[1]) / 2

    ProbePeaks = np.delete(ProbePeaks, [0, 1])
    GreyValuePeaks = np.delete(GreyValuePeaks, [0, 1])

    ProbePeaks = np.insert(ProbePeaks, [0], Corr_peak)
    GreyValuePeaks = np.insert(GreyValuePeaks, [0], Corr_Grey_Value)

    # Get Distance from peaks
    ProbePeaks = ProbePeaks * (18 / (len(Probe.Pixel) + 6))
    LightPeaks = LightPeaks * (18 / (len(Probe.Pixel) + 6))

    # GreyValuePeaks = GreyValuePeaks[:11]
    # ProbePeaks = ProbePeaks[:11]

    # Plot peaks
    plt.figure()
    plt.plot(Probe.Distance, Probe.GreyValue, ls='--', color='blue',
             label="Grauwert")
    plt.plot(ProbePeaks, GreyValuePeaks, color='black',
             ls='', marker='o', label="Erkannte Peaks")
    # plt.plot(LightPeaks, GreyValueLightPeaks, color='grey',
    #          ls='', marker='o', label="Erkannte, schwache Peaks")
    plt.xlabel(r"$r / \mathrm{cm}$")
    plt.ylabel('inverser Grauwert')
    plt.xlim(0, 18)
    plt.legend(loc="lower left")
    plt.tight_layout
    plt.savefig("Auswertung/Grafiken/" + Probe.name + "_Peaks.pdf")

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
    print(PeakAngle)
    PeakAngle = np.sort(PeakAngle)
    print(PeakAngle)

    # convert the angles to interplanar distance according to braggs law
    d = lam / (2 * np.sin(0.5 * PeakAngle * np.pi / 180))

    print("ZnS n=sqrt(h**2+k**2+l**2):")
    ZnS_n, ZnS_a = find_lattice_constants(d, 'ZnS', 7)
    print(ZnS_n)
    print('ZnS h, k, l:')
    ZnS_h, ZnS_k, ZnS_l = find_hkl('ZnS', ZnS_n, 7)
    print(ZnS_h, ZnS_k, ZnS_l)
    print("ZnS a:")
    print(ZnS_a)
    print("CsCl n=sqrt(h**2+k**2+l**2):")
    CsCl_n, CsCl_a = find_lattice_constants(d, 'CsCl', 7)
    print(CsCl_n)
    print('CsCl h, k, l:')
    CsCl_h, CsCl_k, CsCl_l = find_hkl('CsCl', CsCl_n, 7)
    print(CsCl_h, CsCl_k, CsCl_l)
    print("CsCl a:")
    print(CsCl_a)
    print("NaCl n=sqrt(h**2+k**2+l**2):")
    NaCl_n, NaCl_a = find_lattice_constants(d, 'NaCl', 7)
    print(NaCl_n)
    print('NaCl h, k, l:')
    NaCl_h, NaCl_k, NaCl_l = find_hkl('NaCl', NaCl_n, 7)
    print(NaCl_h, NaCl_k, NaCl_l)
    print("NaCl a:")
    print(NaCl_a)
    print("F n=sqrt(h**2+k**2+l**2):")
    F_n, F_a = find_lattice_constants(d, 'F', 7)
    print(F_n)
    print('F h, k, l:')
    F_h, F_k, F_l = find_hkl('F', F_n, 7)
    print(F_h, F_k, F_l)
    print("F a:")
    print(F_a)

    np.savetxt("Auswertung/Grafiken/" + Probe.name +
               "_ZnS_CsCl_Tabelle.tex",
               np.column_stack([
                               PeakAngle,
                               ZnS_n**2,
                               ZnS_h,
                               ZnS_k,
                               ZnS_l,
                               ZnS_a * 10**(12),
                               CsCl_n**2,
                               CsCl_h,
                               CsCl_k,
                               CsCl_l,
                               CsCl_a * 10**(12),
                               ]), delimiter=' & ', newline=r' \\' + '\n',
               fmt='%.2f & %.0f & %.0f & %.0f & %.0f & %.2f & %.0f & %.0f & %.0f & %.0f & %.2f')
    np.savetxt("Auswertung/Grafiken/" + Probe.name +
               "_NaCl_F_Tabelle.tex",
               np.column_stack([
                               PeakAngle,
                               NaCl_n**2,
                               NaCl_h,
                               NaCl_k,
                               NaCl_l,
                               NaCl_a * 10**(12),
                               F_n**2,
                               F_h,
                               F_k,
                               F_l,
                               F_a * 10**(12),
                               ]), delimiter=' & ', newline=r' \\' + '\n',
               fmt='%.2f & %.0f & %.0f & %.0f & %.0f & %.2f & %.0f & %.0f & %.0f & %.0f & %.2f')

    # Compute best_fit for a thrugh linear regression
    ZnS_params, cov = curve_fit(linear, np.cos(0.5 * PeakAngle * np.pi /
                                180)**2, ZnS_a)
    ZnS_errors = np.sqrt(np.diag(cov))
    m_ZnS = ufloat(ZnS_params[0], ZnS_errors[0])
    n_ZnS = ufloat(ZnS_params[1], ZnS_errors[1])

    CsCl_params, cov = curve_fit(linear, np.cos(0.5 * PeakAngle * np.pi /
                                180)**2, CsCl_a)
    CsCl_errors = np.sqrt(np.diag(cov))
    m_CsCl = ufloat(CsCl_params[0], CsCl_errors[0])
    n_CsCl = ufloat(CsCl_params[1], CsCl_errors[1])

    NaCl_params, cov = curve_fit(linear, np.cos(0.5 * PeakAngle * np.pi /
                                180)**2, NaCl_a)
    NaCl_errors = np.sqrt(np.diag(cov))
    m_NaCl = ufloat(NaCl_params[0], NaCl_errors[0])
    n_NaCl = ufloat(NaCl_params[1], NaCl_errors[1])

    F_params, cov = curve_fit(linear, np.cos(0.5 * PeakAngle * np.pi /
                                180)**2, F_a)
    F_errors = np.sqrt(np.diag(cov))
    m_F = ufloat(F_params[0], F_errors[0])
    n_F = ufloat(F_params[1], F_errors[1])

    print("ZnS best fit:")
    print("m = ", m_ZnS, ", n = ", n_ZnS)
    print("CsCl best fit:")
    print("m = ", m_CsCl, ", n = ", n_CsCl)
    print("NaCl best fit:")
    print("m = ", m_NaCl, ", n = ", n_NaCl)
    print("F best fit:")
    print("m = ", m_F, ", n = ", n_F)

    x_range = np.linspace(0, 90, 1000)
    x_range = np.cos(x_range * np.pi / 180)**2

    # Plot peaks
    plt.figure()
    plt.plot(np.cos(PeakAngle * 0.5 * np.pi / 180)**2, ZnS_a * 10**(12),
             marker='x', color='blue', ls='')
    plt.plot(x_range, linear(x_range, *ZnS_params) * 10**(12),
             ls='-', color='blue', label='Hypothese: ZnS-Gitter')
    plt.xlabel(r"$\cos{(\phi)}^{2}$")
    plt.ylabel(r'Berechnete Gitterkonstante$ / \mathrm{pm}$')
    plt.legend(loc="best")
    plt.xlim(0, 1)
    plt.tight_layout
    plt.savefig("Auswertung/Grafiken/" +
                Probe.name + "_ZnS_Ausgleichsrechnung.pdf")

    plt.figure()
    plt.plot(np.cos(PeakAngle * 0.5 * np.pi / 180)**2, CsCl_a * 10**(12),
             marker='x', color='green', ls='')
    plt.plot(x_range, linear(x_range, *CsCl_params) * 10**(12),
             ls='-', color='green', label='Hypothese: CsCl-Gitter')
    plt.xlabel(r"$\cos{(\phi)}^{2}$")
    plt.ylabel(r'Berechnete Gitterkonstante$ / \mathrm{pm}$')
    plt.legend(loc="best")
    plt.xlim(0, 1)
    plt.tight_layout
    plt.savefig("Auswertung/Grafiken/" +
                Probe.name + "_CsCl_Ausgleichsrechnung.pdf")

    plt.figure()
    plt.plot(np.cos(PeakAngle * 0.5 * np.pi / 180)**2, NaCl_a * 10**(12),
             marker='x', color='orange', ls='')
    plt.plot(x_range, linear(x_range, *NaCl_params) * 10**(12),
             ls='-', color='orange', label='Hypothese: NaCl-Gitter')
    plt.xlabel(r"$\cos{(\phi)}^{2}$")
    plt.ylabel(r'Berechnete Gitterkonstante$ / \mathrm{pm}$')
    plt.legend(loc="best")
    plt.xlim(0, 1)
    plt.tight_layout
    plt.savefig("Auswertung/Grafiken/" +
                Probe.name + "_NaCl_Ausgleichsrechnung.pdf")

    plt.figure()
    plt.plot(np.cos(PeakAngle * 0.5 * np.pi / 180)**2, F_a * 10**(12),
             marker='x', color='black', ls='')
    plt.plot(x_range, linear(x_range, *F_params) * 10**(12),
             ls='-', color='black', label='Hypothese: F-Gitter')
    plt.xlabel(r"$\cos{(\phi)}^{2}$")
    plt.ylabel(r'Berechnete Gitterkonstante$ / \mathrm{pm}$')
    plt.legend(loc="best")
    plt.xlim(0, 1)
    plt.tight_layout
    plt.savefig("Auswertung/Grafiken/" +
                Probe.name + "_F_Ausgleichsrechnung.pdf")

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(np.cos(PeakAngle * 0.5 * np.pi / 180)**2, ZnS_a * 10**(12),
            marker='x', color='blue', ls='')
    ax.plot(x_range, linear(x_range, *ZnS_params) * 10**(12),
            ls='-', color='blue', label='Hypothese: ZnS-Gitter')
    ax.plot(np.cos(PeakAngle * 0.5 * np.pi / 180)**2, CsCl_a * 10**(12),
            marker='x', color='green', ls='')
    ax.plot(x_range, linear(x_range, *CsCl_params) * 10**(12),
            ls='-', color='green', label='Hypothese: CsCl-Gitter')
    ax.plot(np.cos(PeakAngle * 0.5 * np.pi / 180)**2, NaCl_a * 10**(12),
            marker='x', color='orange', ls='')
    ax.plot(x_range, linear(x_range, *NaCl_params) * 10**(12),
            ls='-', color='orange', label='Hypothese: NaCl-Gitter')
    ax.plot(np.cos(PeakAngle * 0.5 * np.pi / 180)**2, F_a * 10**(12),
            marker='x', color='black', ls='')
    ax.plot(x_range, linear(x_range, *F_params) * 10**(12),
            ls='-', color='black', label='Hypothese: F-Gitter')
    ax.set_xlabel(r"$\cos{(\phi)}^{2}$")
    ax.set_ylabel(r'Berechnete Gitterkonstante$ / \mathrm{pm}$')
    art = []
    lgd = pylab.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2)
    art.append(lgd)
    ax.set_xlim(0, 1)
    plt.tight_layout
    pylab.savefig("Auswertung/Grafiken/" +
                  Probe.name + "_Ausgleichsrechnung.pdf",
                  additional_artists=art,
                  bbox_inches="tight")

print("------------------------------------------------------------------")
print('Thats all folks!')
