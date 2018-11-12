import numpy as np

# Possible Millerindizes for fcc and bcc
############
# 1: bcc: h+k+l even
# 2: fcc: h,k,l all even or odd
############
# Therefor use the the modulo operation: n % 2 == {0: n even or 1: n odd}

# Set maximum-value for witch the Millerindizes should be computed, meaning the
# maximal value for h, k, and l

max_value = 7

# initialize lists fpr indizes
h_bcc = []
k_bcc = []
l_bcc = []

h_fcc = []
k_fcc = []
l_fcc = []

# bcc
print("bcc:")
for h in range(max_value):
	for k in range(max_value):
		for l in range(max_value):
			if (h+k+l) % 2 == 0 and (h+k+l) > 0:
				print(h, k, l)
				h_bcc.append(h)
				k_bcc.append(k)
				l_bcc.append(l)

# fcc
print("fcc:")
for h in range(max_value):
	for k in range(max_value):
		for l in range(max_value):
			if (h % 2 + k % 2 + l % 2) == 0 or (h % 2 + k % 2 + l % 2) == 1:
				if (h+k+l) > 0:
					print(h, k, l)
					h_fcc.append(h)
					k_fcc.append(k)
					l_fcc.append(l)
