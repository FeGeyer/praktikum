import numpy as np
import matplotlib.pyplot as plt
import millerindizes as miller

h, k, l = miller.CsCl(10)

x_Dia = np.array([0, 0.5, 0.5, 0  , 0.25, 0.75, 0.75, 0.25])
y_Dia = np.array([0, 0.5, 0  , 0.5, 0.25, 0.75, 0.25, 0.75])
z_Dia = np.array([0, 0  , 0.5, 0.5, 0.25, 0.25, 0.75, 0.75])

x_NaCl = np.array([0, 0.5, 0.5, 0  , 0.5, 1  , 1   , 0.5])
y_NaCl = np.array([0, 0.5, 0  , 0.5, 0.5, 1  , 0.25, 1  ])
z_NaCl = np.array([0, 0  , 0.5, 0.5, 0.5, 0.5, 1   , 1  ])

x_CsCl = np.array([0, 0.5])
y_CsCl = np.array([0, 0.5])
z_CsCl = np.array([0, 0.5])

x_F = np.array([0, 0.5, 0.5, 0  , 0.25, 0.75, 0.75, 0.25, 0.75, 0.25, 0.25, 0.75])
y_F = np.array([0, 0.5, 0  , 0.5, 0.25, 0.75, 0.25, 0.75, 0.75, 0.25, 0.75, 0.25])
z_F = np.array([0, 0  , 0.5, 0.5, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75, 0.25, 0.25])


def S(h, k, l, x, y, z):
    return np.exp((-1) * np.pi * 1j * (h * x + k * y + l * z))


def Sbar(h, k, l, x, y, z):
    return np.exp((1) * np.pi * 1j * (h * x + k * y + l * z))

n = h**2 + k**2 + l**2
S2_ZnS = np.empty([1])
S2_NaCl = np.empty([1])
S2_CsCl = np.empty([1])
S2_F = np.empty([1])

for i in range(len(h)):
    print(h[i], k[i], l[i])
    S2_ZnS = np.append(S2_ZnS, np.real(np.sum(S(h[i], k[i], l[i], x_Dia, y_Dia, z_Dia)) *
             np.sum(Sbar(h[i], k[i], l[i], x_Dia, y_Dia, z_Dia))))
    S2_NaCl = np.append(S2_NaCl, np.real(np.sum(S(h[i], k[i], l[i], x_NaCl, y_NaCl, z_NaCl)) *
              np.sum(Sbar(h[i], k[i], l[i], x_NaCl, y_NaCl, z_NaCl))))
    S2_CsCl = np.append(S2_CsCl, np.real(np.sum(S(h[i], k[i], l[i], x_CsCl, y_CsCl, z_CsCl)) *
              np.sum(Sbar(h[i], k[i], l[i], x_CsCl, y_CsCl, z_CsCl))))
    S2_F = np.append(S2_F, np.real(np.sum(S(h[i], k[i], l[i], x_F, y_F, z_F)) *
           np.sum(Sbar(h[i], k[i], l[i], x_F, y_F, z_F))))

S2_ZnS = np.delete(S2_ZnS, [0])
S2_NaCl = np.delete(S2_NaCl, [0])
S2_CsCl = np.delete(S2_CsCl, [0])
S2_F = np.delete(S2_F, [0])


# Plot peaks
plt.figure()
plt.plot(n, S2_ZnS, marker='x', color='red', ls='-')
plt.xlabel(r"$h^{2} + k^{2} + l^{2}$")
plt.ylabel(r'$|S|^{2}$')
plt.legend(loc="best")
plt.tight_layout
plt.show()

plt.figure()
plt.plot(n, S2_NaCl, marker='x', color='red', ls='-')
plt.xlabel(r"$h^{2} + k^{2} + l^{2}$")
plt.ylabel(r'$|S|^{2}$')
plt.legend(loc="best")
plt.tight_layout
plt.show()

plt.figure()
plt.plot(n, S2_CsCl, marker='x', color='red', ls='-')
plt.xlabel(r"$h^{2} + k^{2} + l^{2}$")
plt.ylabel(r'$|S|^{2}$')
plt.legend(loc="best")
plt.tight_layout
plt.show()

plt.figure()
plt.plot(n, S2_F, marker='x', color='red', ls='-')
plt.xlabel(r"$h^{2} + k^{2} + l^{2}$")
plt.ylabel(r'$|S|^{2}$')
plt.legend(loc="best")
plt.tight_layout
plt.show()
