import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#
# LEGO plots {{{

#  TCAG
# T
# C
# A
# G
# 
# C>G C>A C>T   4 3 5
# A>T A>C A>G   2 0 1

# sort ch96 by row-major ordering according to LEGO layout
y = np.r_[0:96];

L = pd.DataFrame({ "ch" : np.right_shift(y & 0x70, 4), "l" : np.right_shift(y & 0xC, 2), "r" : y & 3})
L["ch_c"] = pd.Categorical(L["ch"], [4, 3, 5, 2, 0, 1]);
L["l_c"] = pd.Categorical(L["l"], [3, 1, 0, 2]);
L["r_c"] = pd.Categorical(L["r"], [3, 1, 0, 2]);

ch96_grid = np.hstack([
              np.reshape(L.loc[L["ch"] > 2, :] \
                         .sort_values(["ch_c", "r_c", "l_c"]).index.values, (-1, 4)),
              np.reshape(L.loc[L["ch"] < 3, :] \
                         .sort_values(["ch_c", "r_c", "l_c"]).index.values, (-1, 4))
]).T

lego_colors = np.array([0,   0.2, 0.8,
                        0.1, 0.8, 0.1,
                        0.5, 0.3, 0.7, 
                        0,   0.7, 0.7,
                        1,   0,   0,
                        1,   1,   0]).reshape(-1, 3);
lego_colors = lego_colors[np.right_shift(ch96_grid & 0x70, 4)];

def lego(ch96_counts, fnum = None, axes = None):
	xc, yc = np.meshgrid(np.r_[0:12], np.r_[0:8]);

	if axes is None:
		if fnum is None:
			f = plt.figure(); plt.clf()
			fnum = f.number
		else:
			f = plt.figure(fnum)

		ax = f.add_subplot(projection = '3d', proj_type = 'ortho', azim = 60)
	else:
		ax = axes

	# TODO: specify z-index of each bar individually
	#np.fliplr(xc + np.r_[0:88:11][:, None] + yc)

	ax.bar3d(xc.ravel(), yc.ravel(), 0, 0.8, 0.8, ch96_counts[ch96_grid].ravel(), lego_colors.reshape([96, 3, -1]).squeeze(), edgecolor = 'k', shade = False, zsort = 'max')
	ax.set_xlim3d(12.4, -0.6)
	ax.set_xticks([])
	ax.set_yticks([])

	return ax

# }}}

def logticks(mn, mx):
	minb = np.floor(np.log10(mn))
	maxb = np.floor(np.log10(mx))

	t = [];
	for i in np.r_[minb:(maxb + 1)]:
		t.extend(np.r_[1:10]*10.0**i)

	t = np.array(t); t = t[~((t < mn) | (t > mx))];

	return t

#
# spine manipulation {{{

def hide_spines(ax = None, offset = None):
	if ax is None:
		ax = plt.gca()
	if offset is None:
		offset = -0.03

	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)

	ax.spines["left"].set_position(("axes", offset))
	ax.spines["bottom"].set_position(("axes", offset))

def spine_bounds(ax = None, t = None, r = None, b = None, l = None): 
	if ax is None:
		ax = plt.gca()

	for s, n in zip([t, r, b, l], ["top", "right", "bottom", "left"]):
		if s is None:
			continue

		ax.spines[n].set_bounds(*s)

# }}}
