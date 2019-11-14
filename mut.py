from pyfaidx import Fasta as _Fasta
import numpy as _np
import sys as _sys
import mmap as _mmap
import pandas as _pd
import tqdm as _tqdm
from . import seq as _seq

_byte_LuT = [None]*256
for i in _np.r_[0:256].astype(_np.uint8):
	_byte_LuT[i] = _np.flatnonzero(_np.unpackbits(i))

def filter_mutations_against_gnomAD(M, ref = None, field_map = None, gnomad_dir = None):
	# set gnomAD bitwise track reference directory
	_seq.set_gnomad_ref_params(gnomad_dir = gnomad_dir)

	# assume default columns "chr", "pos", "ref", "newbase"
	if field_map is None:
		field_map = dict(zip(M.columns, range(0, len(M.columns))))

	if not _np.all([x in field_map for x in ["chr", "pos", "ref", "newbase"]]):
		raise KeyError("Default fieldname not found in columns!")

	# get column names of chromosome, position, reference/mutant alleles
	chr_f = M.columns[field_map["chr"]]
	pos_f = M.columns[field_map["pos"]]
	ref_f = M.columns[field_map["ref"]]
	newbase_f = M.columns[field_map["newbase"]]

	if not _np.issubdtype(M[chr_f].dtype, _np.integer):
		raise ValueError("Chromosomes must be specified as integers, 1-24")

	M["SSNV_idx"] = M[ref_f].isin(list("ACGT")) & M[newbase_f].isin(list("ACGT"))

	# store overlapping coordinates in dataframe
	O = []
	
	# loop over chromosomes; query gnomAD for each
	for ch, Mc in _tqdm.tqdm(M.groupby(chr_f)):
		# gnomAD has no sex chromosome variants
		if ch >= 23:
			break

		z = _np.zeros(_seq.get_chrlens(ref = ref)[ch - 1], dtype = _np.bool);

		# get 1 bit packed gnomAD representations for:
		# - all alleles
		# - SNP-specific alleles
		for base in ["all"] + list("ACGT"):
			m = _np.copy(z)
			if base != "all":
				_seq.set_gnomad_ref_params(bin_stem = "to_" + base)

				m[Mc.loc[Mc["SSNV_idx"] & (Mc[newbase_f] == base), pos_f] - 1] = True;
			else:
				m[Mc[pos_f] - 1] = True;

			# pack mutation coordinates into 1 bit array
			m = _np.packbits(m)

			# get gnomAD packed array
			g = _seq.query_gnomad_1bit_raw(ch);

			# compute bitwise intersection; index nonzero sites
			bitwise_overlap = g & m 
			bwol_idx = _np.flatnonzero(bitwise_overlap)

			# unpack bitwise intersection to coordinates and append to DF
			ol_pos = _np.hstack([
					_byte_LuT[byte] + 8*idx for byte, idx in
			        zip(bitwise_overlap[bwol_idx], bwol_idx)
			]) + 1
			O.append(_pd.DataFrame({ "chr" : ch, "pos" : ol_pos, "allele" : base }))

	O = _pd.concat(O)

	# intersect with mutations
	M["gpos"] = _seq.chrpos2gpos(M["chr"], M["pos"], ref = ref)
	O["gpos"] = _seq.chrpos2gpos(O["chr"], O["pos"], ref = ref)

	for base in ["all"] + list("ACGT"):
		M["gnomAD_" + base] = M["gpos"].isin(O.loc[O["allele"] == base, "gpos"])

	return M.drop(columns = ["gpos", "SSNV_idx"])