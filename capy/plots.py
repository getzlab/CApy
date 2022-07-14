import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from statsmodels.stats import multitest
from adjustText import adjust_text

from capy import num as num

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

	ax.view_init(54, -124) 
	for i, (x, y, z, color) in enumerate(zip(
	  np.fliplr(xc).ravel(),
	  np.fliplr(np.flipud(yc)).ravel(),
	  np.fliplr(ch96_counts[ch96_grid]).ravel(),
	  np.fliplr(lego_colors).reshape([96, 3, -1]).squeeze()
	)):
		b3 = ax.bar3d(x, y, 0, 0.8, 0.8, z, color, edgecolor = 'k', shade = False, zsort = "max")
		b3.set_sort_zpos(i)
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

#
# Q-Q

def QQ(pvalues, labels = None, sig_thresh = 0.1, near_sig_thresh = 0.25, fnum = None):
    si = np.argsort(-np.log10(pvalues))
    logp = -np.log10(pvalues)[si]
    n_p = len(logp)

    # expected quantiles
    x = -np.log10(np.r_[n_p:0:-1]/(n_p + 1))

    #
    # FDR
    _, q, _, _ = multitest.multipletests(pvalues[si], method = "fdr_bh")
    sig_idx = q < sig_thresh
    n_sig_idx = q < near_sig_thresh

    fdr_color_idx = np.c_[sig_idx, n_sig_idx]@np.r_[2, 1]

    #                      nonsig     near sig   <padding>     sig
    fdr_colors = np.array([[0, 0, 1], [0, 1, 1], [-1, -1, -1], [1, 0, 0]])

    #
    # plot
    f = plt.figure(fnum); plt.clf()
    plt.scatter(x, logp, c = fdr_colors[fdr_color_idx])

    # 1:1 line
    plt.plot(plt.xlim(), plt.xlim(), color = 'k', linestyle = ':')

    # confidence interval (TODO)

    #
    # labels (if given)
    if labels is not None:
        if len(labels) != len(x):
            raise ValueError("Length of labels must match length of p-value array")
        labels = labels[si]
        label_plt = [plt.text(x, y, l) for x, y, l in zip(x[n_sig_idx], logp[n_sig_idx], labels[n_sig_idx])]
        adjust_text(label_plt, arrowprops = { 'color' : 'k', "arrowstyle" : "-" })

    plt.xlabel("Expected quantile (-log10)")
    plt.ylabel("Observed quantile (-log10)")

    return f

#
# pixel scatterplot
def pixplot(x, y, **kwargs):
    return plt.plot(x, y, marker = ',', linewidth = 0, **kwargs)
