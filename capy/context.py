import numpy as np
import pandas as pd

B = dict(zip(list("ACGT"), range(0, 4)));

# lookup table for strand-collapsing c65
x = (~np.r_[0:32]) & 0x3F;
x = (x & 0x30) | np.left_shift(x & 0x3, 2) | np.right_shift(x & 0xC, 2)
c32LuT = np.full(64, 0, dtype = np.uint8);
c32LuT[x] = np.r_[0:32];

def c65_to_ch96(c65, newbase):
	c65 = np.array(c65);
	nbidx = np.array([B[x] for x in newbase]);

	# strand collapse GT -> CA
	GT_idx = (c65 & 0x20) > 0;
	c32 = c65;
	c32[GT_idx] = c32LuT[c32[GT_idx]]; 
	nbidx[GT_idx] = 3 - nbidx[GT_idx];

	# convert newbases {A->[CGT] / C->[ACT]} -> [012]
	A_idx = c32 < 16;
	nbidx[A_idx] = nbidx[A_idx] - 1;
	nbidx[~A_idx & (nbidx >= 2)] = nbidx[~A_idx & (nbidx >= 2)] - 1;

	return ((~A_idx)*48 + 16*nbidx) | (c32 & 0xF);

def c65_to_c32(c65):
	GT_idx = (c65 & 0x20) > 0;
	c32 = c65;
	c32[GT_idx] = c32LuT[c32[GT_idx]]; 
	return c32

# strings for printing context
bc = np.array(['A->C', 'A->G', 'A->T', 'C->A', 'C->G', 'C->T']);
b = np.array(list("ACGT"));
y = np.r_[0:96];
ch96_def = pd.Series(["".join(x) for x in
                      zip(b[np.right_shift(y & 0xC, 2)], 96*["("],
                          bc[np.right_shift(y & 0x70, 4)], 96*[")"], b[y & 3])])

def print_ch96():
	return ch96_def
