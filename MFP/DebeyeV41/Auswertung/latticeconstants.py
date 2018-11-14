import numpy as np
# Possible Millerindizes for fcc and bcc
############
# 1: bcc: h+k+l even
# 2: fcc: h,k,l all even or odd
############
# Therefor use the the modulo operation: n % 2 == {0: n even or 1: n odd}
# max_value denotes the maximal value for h, k and l respectivly that should
# be evaluated


def find_permutations(array):
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
	array = np.empty([len(h), 3])
	for i in range(len(h)):
		row = -np.array([h[i], k[i], l[i]])
		# print(np.sort(row))
		array[i, :] = -np.sort(row)
	return array


# bcc
def bcc(max_value):
	h_bcc = np.empty([1])
	k_bcc = np.empty([1])
	l_bcc = np.empty([1])
	print("bcc:")
	for l in range(max_value):
		for k in range(max_value):
			for h in range(max_value):
				if (h + k + l) % 2 == 0 and (h + k + l) > 0:
					h_bcc = np.append(h_bcc, h)
					k_bcc = np.append(k_bcc, k)
					l_bcc = np.append(l_bcc, l)
	h_bcc = np.delete(h_bcc, [0])
	k_bcc = np.delete(k_bcc, [0])
	l_bcc = np.delete(l_bcc, [0])

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

	indizes = find_permutations(combined_array[:, 0:3])

	h_bcc = np.delete(combined_array[:, 0], indizes)
	k_bcc = np.delete(combined_array[:, 1], indizes)
	l_bcc = np.delete(combined_array[:, 2], indizes)

	combined_array = np.empty([len(h_bcc), 4])
	combined_array[:, 0] = h_bcc
	combined_array[:, 1] = k_bcc
	combined_array[:, 2] = l_bcc
	combined_array[:, 3] = np.sqrt(h_bcc**2 + k_bcc**2 + l_bcc**2)

	combined_array = combined_array[combined_array[:, 3].argsort()]

	return combined_array[:, 0], combined_array[:, 1], combined_array[:, 2]


# fcc
def fcc(max_value):
	h_fcc = np.empty([])
	k_fcc = np.empty([])
	l_fcc = np.empty([])
	print("fcc:")
	for l in range(max_value):
		for k in range(max_value):
			for h in range(max_value):
				if (h % 2 + k % 2 + l % 2) == 0 or (h % 2 + k % 2 + l % 2) == 1:
					if (h + k + l) > 0:
						h_fcc = np.append(h_fcc, h)
						k_fcc = np.append(k_fcc, k)
						l_fcc = np.append(l_fcc, l)
	h_fcc = np.delete(h_fcc, [0])
	k_fcc = np.delete(k_fcc, [0])
	l_fcc = np.delete(l_fcc, [0])

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

	print(combined_array)

	indizes = find_permutations(combined_array[:, 0:3])

	h_fcc = np.delete(combined_array[:, 0], indizes)
	k_fcc = np.delete(combined_array[:, 1], indizes)
	l_fcc = np.delete(combined_array[:, 2], indizes)

	combined_array = np.empty([len(h_fcc), 4])
	combined_array[:, 0] = h_fcc
	combined_array[:, 1] = k_fcc
	combined_array[:, 2] = l_fcc
	combined_array[:, 3] = np.sqrt(h_fcc**2 + k_fcc**2 + l_fcc**2)

	combined_array = combined_array[combined_array[:, 3].argsort()]

	print(combined_array)

	return combined_array[:, 0], combined_array[:, 1], combined_array[:, 2]
