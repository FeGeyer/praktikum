import numpy as np
'''
Some functions to compute possible miller-tupel (h, k, l) for fcc-
and bcc-lattice.

Fuctions make use of the following selection rules:
1: bcc: h+k+l even
2: fcc: h,k,l all even or odd
3: Diamant: (h+k+l)/2 even or h+k+l odd
4: CsCl: all h,k,l allowed
Therefor use the the modulo operation: n % 2 == {0: n even or 1: n odd}.
'''


def find_permutations(array):

	'''
	3darray: array of ints, array of the sorted millerindizes h, k, l

	Funktion to find the index of permutated or doubled miller-tupel (h, k, l).
	Needs row-wise sorted (h, k, l).

	Returns the indizes of permutations and doubles.
	'''

	indize = np.empty([1])
	for i in range(len(array)):
		tested = array[i, :]
		for j in range(i, len(array) - 1):
			next = array[j + 1, :]
			if tested[0] == next[0] and tested[1] == next[1] and tested[2] == next[2]:
				indize = np.append(indize, j + 1)
	indize = np.delete(indize, [0])
	return indize


def sort_rows(h, k, l):

	'''
	h, k, l: arrays of ints, arrays millerindizes.

	Funktion to sort miller-tupels (h, k, l) row-wise.

	Returns 3darray of row-wise sorted miller-tupel.
	'''

	array = np.empty([len(h), 3])
	for i in range(len(h)):
		row = -np.array([h[i], k[i], l[i]])
		# print(np.sort(row))
		array[i, :] = -np.sort(row)
	return array


# bcc
def bcc(max_value):

	'''
	max_value: Maximal value for h, k, l respectively.

	Funktion compute possible bcc millerindizes. Calls functions to sort out
	permutations and doubles.

	Returns arrays of ints with millerindizes for the bcc-lattice.
	'''

	# initialize arrays
	h_bcc = np.empty([1])
	k_bcc = np.empty([1])
	l_bcc = np.empty([1])

	# Compute all h, k, l, witch satisfy the primary criterion for bcc-lattice
	for l in range(max_value):
		for k in range(max_value):
			for h in range(max_value):
				if (h + k + l) % 2 == 0 and (h + k + l) > 0:
					h_bcc = np.append(h_bcc, h)
					k_bcc = np.append(k_bcc, k)
					l_bcc = np.append(l_bcc, l)

	# delete first initialized value
	h_bcc = np.delete(h_bcc, [0])
	k_bcc = np.delete(k_bcc, [0])
	l_bcc = np.delete(l_bcc, [0])

	# sort
	array = sort_rows(h_bcc, k_bcc, l_bcc)
	h_bcc = array[:, 0]
	k_bcc = array[:, 1]
	l_bcc = array[:, 2]

	combined_array = np.empty([len(h_bcc), 4])
	combined_array[:, 0] = h_bcc
	combined_array[:, 1] = k_bcc
	combined_array[:, 2] = l_bcc
	combined_array[:, 3] = np.sqrt(h_bcc**2 + k_bcc**2 + l_bcc**2)

	combined_array = combined_array[combined_array[:, 3].argsort()]

	# find permutations an delete them
	indizes = find_permutations(combined_array[:, 0:3])

	h_bcc = np.delete(combined_array[:, 0], indizes)
	k_bcc = np.delete(combined_array[:, 1], indizes)
	l_bcc = np.delete(combined_array[:, 2], indizes)

	combined_array = np.empty([len(h_bcc), 4])
	combined_array[:, 0] = h_bcc
	combined_array[:, 1] = k_bcc
	combined_array[:, 2] = l_bcc
	combined_array[:, 3] = np.sqrt(h_bcc**2 + k_bcc**2 + l_bcc**2)

	# sort
	combined_array = combined_array[combined_array[:, 3].argsort()]

	return combined_array[:, 0], combined_array[:, 1], combined_array[:, 2]


# fcc
def fcc(max_value):

	'''
	max_value: Maximal value for h, k, l respectively.

	Funktion compute possible fcc millerindizes. Calls functions to sort out
	permutations and doubles.

	Returns arrays of ints with millerindizes for the fcc-lattice.
	'''

	# initialize arrays
	h_fcc = np.empty([])
	k_fcc = np.empty([])
	l_fcc = np.empty([])

	# Compute all h, k, l, witch satisfy the primary criterion for fcc-lattice
	for l in range(max_value):
		for k in range(max_value):
			for h in range(max_value):
				if ((h % 2) + (k % 2) + (l % 2)) == 0 or ((h % 2) + (k % 2) + (l % 2)) == 3:
					if (h + k + l) > 0:
						h_fcc = np.append(h_fcc, h)
						k_fcc = np.append(k_fcc, k)
						l_fcc = np.append(l_fcc, l)

	# delete first initialized value
	h_fcc = np.delete(h_fcc, [0])
	k_fcc = np.delete(k_fcc, [0])
	l_fcc = np.delete(l_fcc, [0])

	# sort
	array = sort_rows(h_fcc, k_fcc, l_fcc)
	h_fcc = array[:, 0]
	k_fcc = array[:, 1]
	l_fcc = array[:, 2]

	combined_array = np.empty([len(h_fcc), 4])
	combined_array[:, 0] = h_fcc
	combined_array[:, 1] = k_fcc
	combined_array[:, 2] = l_fcc
	combined_array[:, 3] = np.sqrt(h_fcc**2 + k_fcc**2 + l_fcc**2)

	combined_array = combined_array[combined_array[:, 3].argsort()]

	# find permutations and delete
	indizes = find_permutations(combined_array[:, 0:3])

	h_fcc = np.delete(combined_array[:, 0], indizes)
	k_fcc = np.delete(combined_array[:, 1], indizes)
	l_fcc = np.delete(combined_array[:, 2], indizes)

	combined_array = np.empty([len(h_fcc), 4])
	combined_array[:, 0] = h_fcc
	combined_array[:, 1] = k_fcc
	combined_array[:, 2] = l_fcc
	combined_array[:, 3] = np.sqrt(h_fcc**2 + k_fcc**2 + l_fcc**2)

	# sort
	combined_array = combined_array[combined_array[:, 3].argsort()]

	return combined_array[:, 0], combined_array[:, 1], combined_array[:, 2]


# Diamant
def Dia(max_value):

	'''
	max_value: Maximal value for h, k, l respectively.

	Funktion compute possible fcc millerindizes. Calls functions to sort out
	permutations and doubles.

	Returns arrays of ints with millerindizes for the Diamant-lattice.
	'''

	# initialize arrays
	h_fcc = np.empty([])
	k_fcc = np.empty([])
	l_fcc = np.empty([])

	# Compute all h, k, l, witch satisfy the primary criterion for fcc-lattice
	for l in range(max_value):
		for k in range(max_value):
			for h in range(max_value):
				if ((h + k + l) / 2) % 2 == 0 or (h + k + l) % 2 == 1:
					if (h + k + l) > 0:
						h_fcc = np.append(h_fcc, h)
						k_fcc = np.append(k_fcc, k)
						l_fcc = np.append(l_fcc, l)

	# delete first initialized value
	h_fcc = np.delete(h_fcc, [0])
	k_fcc = np.delete(k_fcc, [0])
	l_fcc = np.delete(l_fcc, [0])

	# sort
	array = sort_rows(h_fcc, k_fcc, l_fcc)
	h_fcc = array[:, 0]
	k_fcc = array[:, 1]
	l_fcc = array[:, 2]

	combined_array = np.empty([len(h_fcc), 4])
	combined_array[:, 0] = h_fcc
	combined_array[:, 1] = k_fcc
	combined_array[:, 2] = l_fcc
	combined_array[:, 3] = np.sqrt(h_fcc**2 + k_fcc**2 + l_fcc**2)

	combined_array = combined_array[combined_array[:, 3].argsort()]

	# find permutations and delete
	indizes = find_permutations(combined_array[:, 0:3])

	h_fcc = np.delete(combined_array[:, 0], indizes)
	k_fcc = np.delete(combined_array[:, 1], indizes)
	l_fcc = np.delete(combined_array[:, 2], indizes)

	combined_array = np.empty([len(h_fcc), 4])
	combined_array[:, 0] = h_fcc
	combined_array[:, 1] = k_fcc
	combined_array[:, 2] = l_fcc
	combined_array[:, 3] = np.sqrt(h_fcc**2 + k_fcc**2 + l_fcc**2)

	# sort
	combined_array = combined_array[combined_array[:, 3].argsort()]

	return combined_array[:, 0], combined_array[:, 1], combined_array[:, 2]


# Diamant
def All(max_value):

	'''
	max_value: Maximal value for h, k, l respectively.

	Funktion compute possible millerindizes. Calls functions to sort out
	permutations and doubles.

	Returns arrays of ints with millerindizes for the Diamant-lattice.
	'''

	# initialize arrays
	h_fcc = np.empty([])
	k_fcc = np.empty([])
	l_fcc = np.empty([])

	# Compute all h, k, l, witch satisfy the primary criterion for fcc-lattice
	for l in range(max_value):
		for k in range(max_value):
			for h in range(max_value):
				if (h + k + l) > 0:
					h_fcc = np.append(h_fcc, h)
					k_fcc = np.append(k_fcc, k)
					l_fcc = np.append(l_fcc, l)

	# delete first initialized value
	h_fcc = np.delete(h_fcc, [0])
	k_fcc = np.delete(k_fcc, [0])
	l_fcc = np.delete(l_fcc, [0])

	# sort
	array = sort_rows(h_fcc, k_fcc, l_fcc)
	h_fcc = array[:, 0]
	k_fcc = array[:, 1]
	l_fcc = array[:, 2]

	combined_array = np.empty([len(h_fcc), 4])
	combined_array[:, 0] = h_fcc
	combined_array[:, 1] = k_fcc
	combined_array[:, 2] = l_fcc
	combined_array[:, 3] = np.sqrt(h_fcc**2 + k_fcc**2 + l_fcc**2)

	combined_array = combined_array[combined_array[:, 3].argsort()]

	# find permutations and delete
	indizes = find_permutations(combined_array[:, 0:3])

	h_fcc = np.delete(combined_array[:, 0], indizes)
	k_fcc = np.delete(combined_array[:, 1], indizes)
	l_fcc = np.delete(combined_array[:, 2], indizes)

	combined_array = np.empty([len(h_fcc), 4])
	combined_array[:, 0] = h_fcc
	combined_array[:, 1] = k_fcc
	combined_array[:, 2] = l_fcc
	combined_array[:, 3] = np.sqrt(h_fcc**2 + k_fcc**2 + l_fcc**2)

	# sort
	combined_array = combined_array[combined_array[:, 3].argsort()]

	return combined_array[:, 0], combined_array[:, 1], combined_array[:, 2]


def S(h, k, l, x, y, z):
    return np.exp((-1) * np.pi * 1j * (h * x + k * y + l * z))


def Sbar(h, k, l, x, y, z):
    return np.exp((1) * np.pi * 1j * (h * x + k * y + l * z))


def ZnS(max_value):
    h, k, l = All(max_value)

    x = np.array([0, 0.5, 0.5, 0  , 0.25, 0.75, 0.75, 0.25])
    y = np.array([0, 0.5, 0  , 0.5, 0.25, 0.75, 0.25, 0.75])
    z = np.array([0, 0  , 0.5, 0.5, 0.25, 0.25, 0.75, 0.75])

    n = h**2 + k**2 + l**2
    S2 = np.empty([1])

    for i in range(len(h)):
        S2 = np.append(S2, np.real(np.sum(S(h[i], k[i], l[i], x, y, z)) *
             np.sum(Sbar(h[i], k[i], l[i], x, y, z))))
    S2 = np.delete(S2, [0])

    # write h, k, l in variables
    h_fcc = h[S2 > 3]
    k_fcc = k[S2 > 3]
    l_fcc = l[S2 > 3]

    # sort
    array = sort_rows(h_fcc, k_fcc, l_fcc)
    h_fcc = array[:, 0]
    k_fcc = array[:, 1]
    l_fcc = array[:, 2]

    combined_array = np.empty([len(h_fcc), 4])
    combined_array[:, 0] = h_fcc
    combined_array[:, 1] = k_fcc
    combined_array[:, 2] = l_fcc
    combined_array[:, 3] = np.sqrt(h_fcc**2 + k_fcc**2 + l_fcc**2)

    combined_array = combined_array[combined_array[:, 3].argsort()]

    # find permutations and delete
    indizes = find_permutations(combined_array[:, 0:3])

    h_fcc = np.delete(combined_array[:, 0], indizes)
    k_fcc = np.delete(combined_array[:, 1], indizes)
    l_fcc = np.delete(combined_array[:, 2], indizes)

    combined_array = np.empty([len(h_fcc), 4])
    combined_array[:, 0] = h_fcc
    combined_array[:, 1] = k_fcc
    combined_array[:, 2] = l_fcc
    combined_array[:, 3] = np.sqrt(h_fcc**2 + k_fcc**2 + l_fcc**2)

    # sort
    combined_array = combined_array[combined_array[:, 3].argsort()]

    return combined_array[:, 0], combined_array[:, 1], combined_array[:, 2]


def NaCl(max_value):
    h, k, l = All(max_value)

    x = np.array([0, 0.5, 0.5, 0  , 0.5, 1  , 1   , 0.5])
    y = np.array([0, 0.5, 0  , 0.5, 0.5, 1  , 0.25, 1  ])
    z = np.array([0, 0  , 0.5, 0.5, 0.5, 0.5, 1   , 1  ])

    n = h**2 + k**2 + l**2
    S2 = np.empty([1])

    for i in range(len(h)):
        S2 = np.append(S2, np.real(np.sum(S(h[i], k[i], l[i], x, y, z)) *
             np.sum(Sbar(h[i], k[i], l[i], x, y, z))))
    S2 = np.delete(S2, [0])

    # write h, k, l in variables
    h_fcc = h[S2 > 3]
    k_fcc = k[S2 > 3]
    l_fcc = l[S2 > 3]

    # sort
    array = sort_rows(h_fcc, k_fcc, l_fcc)
    h_fcc = array[:, 0]
    k_fcc = array[:, 1]
    l_fcc = array[:, 2]

    combined_array = np.empty([len(h_fcc), 4])
    combined_array[:, 0] = h_fcc
    combined_array[:, 1] = k_fcc
    combined_array[:, 2] = l_fcc
    combined_array[:, 3] = np.sqrt(h_fcc**2 + k_fcc**2 + l_fcc**2)

    combined_array = combined_array[combined_array[:, 3].argsort()]

    # find permutations and delete
    indizes = find_permutations(combined_array[:, 0:3])

    h_fcc = np.delete(combined_array[:, 0], indizes)
    k_fcc = np.delete(combined_array[:, 1], indizes)
    l_fcc = np.delete(combined_array[:, 2], indizes)

    combined_array = np.empty([len(h_fcc), 4])
    combined_array[:, 0] = h_fcc
    combined_array[:, 1] = k_fcc
    combined_array[:, 2] = l_fcc
    combined_array[:, 3] = np.sqrt(h_fcc**2 + k_fcc**2 + l_fcc**2)

    # sort
    combined_array = combined_array[combined_array[:, 3].argsort()]

    return combined_array[:, 0], combined_array[:, 1], combined_array[:, 2]


def CsCl(max_value):
    h, k, l = All(max_value)

    x = np.array([0, 0.5])
    y = np.array([0, 0.5])
    z = np.array([0, 0.5])

    n = h**2 + k**2 + l**2
    S2 = np.empty([1])

    for i in range(len(h)):
        S2 = np.append(S2, np.real(np.sum(S(h[i], k[i], l[i], x, y, z)) *
             np.sum(Sbar(h[i], k[i], l[i], x, y, z))))
    S2 = np.delete(S2, [0])

    # write h, k, l in variables
    h_fcc = h[S2 > 3]
    k_fcc = k[S2 > 3]
    l_fcc = l[S2 > 3]

    # sort
    array = sort_rows(h_fcc, k_fcc, l_fcc)
    h_fcc = array[:, 0]
    k_fcc = array[:, 1]
    l_fcc = array[:, 2]

    combined_array = np.empty([len(h_fcc), 4])
    combined_array[:, 0] = h_fcc
    combined_array[:, 1] = k_fcc
    combined_array[:, 2] = l_fcc
    combined_array[:, 3] = np.sqrt(h_fcc**2 + k_fcc**2 + l_fcc**2)

    combined_array = combined_array[combined_array[:, 3].argsort()]

    # find permutations and delete
    indizes = find_permutations(combined_array[:, 0:3])

    h_fcc = np.delete(combined_array[:, 0], indizes)
    k_fcc = np.delete(combined_array[:, 1], indizes)
    l_fcc = np.delete(combined_array[:, 2], indizes)

    combined_array = np.empty([len(h_fcc), 4])
    combined_array[:, 0] = h_fcc
    combined_array[:, 1] = k_fcc
    combined_array[:, 2] = l_fcc
    combined_array[:, 3] = np.sqrt(h_fcc**2 + k_fcc**2 + l_fcc**2)

    # sort
    combined_array = combined_array[combined_array[:, 3].argsort()]

    return combined_array[:, 0], combined_array[:, 1], combined_array[:, 2]


def F(max_value):
    h, k, l = All(max_value)

    x = np.array([0, 0.5, 0.5, 0  , 0.25, 0.75, 0.75, 0.25, 0.75, 0.25, 0.25, 0.75])
    y = np.array([0, 0.5, 0  , 0.5, 0.25, 0.75, 0.25, 0.75, 0.75, 0.25, 0.75, 0.25])
    z = np.array([0, 0  , 0.5, 0.5, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75, 0.25, 0.25])

    n = h**2 + k**2 + l**2
    S2 = np.empty([1])

    for i in range(len(h)):
        S2 = np.append(S2, np.real(np.sum(S(h[i], k[i], l[i], x, y, z)) *
             np.sum(Sbar(h[i], k[i], l[i], x, y, z))))
    S2 = np.delete(S2, [0])

    # write h, k, l in variables
    h_fcc = h[S2 > 3]
    k_fcc = k[S2 > 3]
    l_fcc = l[S2 > 3]

    # sort
    array = sort_rows(h_fcc, k_fcc, l_fcc)
    h_fcc = array[:, 0]
    k_fcc = array[:, 1]
    l_fcc = array[:, 2]

    combined_array = np.empty([len(h_fcc), 4])
    combined_array[:, 0] = h_fcc
    combined_array[:, 1] = k_fcc
    combined_array[:, 2] = l_fcc
    combined_array[:, 3] = np.sqrt(h_fcc**2 + k_fcc**2 + l_fcc**2)

    combined_array = combined_array[combined_array[:, 3].argsort()]

    # find permutations and delete
    indizes = find_permutations(combined_array[:, 0:3])

    h_fcc = np.delete(combined_array[:, 0], indizes)
    k_fcc = np.delete(combined_array[:, 1], indizes)
    l_fcc = np.delete(combined_array[:, 2], indizes)

    combined_array = np.empty([len(h_fcc), 4])
    combined_array[:, 0] = h_fcc
    combined_array[:, 1] = k_fcc
    combined_array[:, 2] = l_fcc
    combined_array[:, 3] = np.sqrt(h_fcc**2 + k_fcc**2 + l_fcc**2)

    # sort
    combined_array = combined_array[combined_array[:, 3].argsort()]

    return combined_array[:, 0], combined_array[:, 1], combined_array[:, 2]
