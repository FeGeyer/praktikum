import importlib.util
spec = importlib.util.spec_from_file_location("millerindizes", "Auswertung/millerindizes.py")
miller = importlib.util.module_from_spec(spec)
spec.loader.exec_module(miller)
import numpy as np


def S(h, k, l, x, y, z):
    return np.exp((-1) * np.pi * 1j * (h * x + k * y + l * z))


def Sbar(h, k, l, x, y, z):
    return np.exp((1) * np.pi * 1j * (h * x + k * y + l * z))


def ZnS(max_value):
    h, k, l = miller.CsCl(max_value)

    x = np.array([0, 0.5, 0.5, 0  , 0.25, 0.75, 0.75, 0.25])
    y = np.array([0, 0.5, 0  , 0.5, 0.25, 0.75, 0.25, 0.75])
    z = np.array([0, 0  , 0.5, 0.5, 0.25, 0.25, 0.75, 0.75])

    n = h**2 + k**2 + l**2
    S2 = np.empty([1])

    for i in range(len(h)):
        S2 = np.append(S2, np.real(np.sum(S(h[i], k[i], l[i], x, y, z)) *
             np.sum(Sbar(h[i], k[i], l[i], x, y, z))))
    S2 = np.delete(S2, [0])

    return n, h[S2 > 2.5], k[S2 > 2.5], l[S2 > 2.5]


def NaCl(max_value):
    h, k, l = miller.CsCl(max_value)

    x = np.array([0, 0.5, 0.5, 0  , 0.5, 1  , 1   , 0.5])
    y = np.array([0, 0.5, 0  , 0.5, 0.5, 1  , 0.25, 1  ])
    z = np.array([0, 0  , 0.5, 0.5, 0.5, 0.5, 1   , 1  ])

    n = h**2 + k**2 + l**2
    S2 = np.empty([1])

    for i in range(len(h)):
        S2 = np.append(S2, np.real(np.sum(S(h[i], k[i], l[i], x, y, z)) *
             np.sum(Sbar(h[i], k[i], l[i], x, y, z))))
    S2 = np.delete(S2, [0])

    return n, h[S2 > 2.5], k[S2 > 2.5], l[S2 > 2.5]


def CsCl(max_value):
    h, k, l = miller.CsCl(max_value)

    x = np.array([0, 0.5])
    y = np.array([0, 0.5])
    z = np.array([0, 0.5])

    n = h**2 + k**2 + l**2
    S2 = np.empty([1])

    for i in range(len(h)):
        S2 = np.append(S2, np.real(np.sum(S(h[i], k[i], l[i], x, y, z)) *
             np.sum(Sbar(h[i], k[i], l[i], x, y, z))))
    S2 = np.delete(S2, [0])

    return n, h[S2 > 2.5], k[S2 > 2.5], l[S2 > 2.5]


def F(max_value):
    h, k, l = miller.CsCl(max_value)

    x = np.array([0, 0.5, 0.5, 0  , 0.25, 0.75, 0.75, 0.25, 0.75, 0.25, 0.25, 0.75])
    y = np.array([0, 0.5, 0  , 0.5, 0.25, 0.75, 0.25, 0.75, 0.75, 0.25, 0.75, 0.25])
    z = np.array([0, 0  , 0.5, 0.5, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75, 0.25, 0.25])

    n = h**2 + k**2 + l**2
    S2 = np.empty([1])

    for i in range(len(h)):
        S2 = np.append(S2, np.real(np.sum(S(h[i], k[i], l[i], x, y, z)) *
             np.sum(Sbar(h[i], k[i], l[i], x, y, z))))
    S2 = np.delete(S2, [0])

    return n, h[S2 > 2.5], k[S2 > 2.5], l[S2 > 2.5]
