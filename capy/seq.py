from pyfaidx import Fasta as _Fasta, FastaNotFoundError as _FastaNotFoundError
import numpy as _np
import sys as _sys
import os as _os
import mmap as _mmap
import pandas as _pd

#
# Reference genome functions
# 

# wrap in a class because we only want to remember which reference we're using,
# so we don't have to instantiate it every time genome_region is called
class _FA:
	def __init__(self, ref_fa_file = "/mnt/j/db/hg19/ref/hs37d5.fa"):
		self._set_reference(ref_fa_file)

	def __auto_reset_reference(f):
		def inner(self, *args, ref = None, **kwargs):
			if ref is None:
				ref = self.ref_fa_file

			if ref != self.ref_fa_file:
				self._set_reference(ref)
				print("NOTE: switched reference to {}".format(ref), file = _sys.stderr)

			return f(self, *args, ref = ref, **kwargs)
		return inner

	@__auto_reset_reference
	def _genome_region(self, chr, start, end, ref = None):
		# TODO: handle implicit chr conversion better
		return self.ref_fa_obj[str(chr)][start:end].seq

	@__auto_reset_reference
	def _chrpos2gpos(self, chr, pos, ref = None):
		# TODO: handle implicit chr conversion better
		#       currently, this only works with b37-style names, and assumes
		#       chr1 comes first
		l = _np.cumsum(_np.r_[0, self.chrlens])
		return l[chr - 1] + pos

	@__auto_reset_reference
	def _gpos2chrpos(self, gpos, ref = None):
		# TODO: handle implicit chr conversion better
		#       currently, this only works with b37-style names, and assumes
		#       chr1 comes first
		l = _np.cumsum(_np.r_[0, self.chrlens])
		chridx = _np.digitize(gpos, l);
		return (chridx, gpos - l[chridx - 1])

	@__auto_reset_reference
	def _get_chrlens(self, ref = None):
		return self.chrlens

	def _set_reference(self, ref):
		try:
			self.ref_fa_file = ref
			self.ref_fa_obj = _Fasta(self.ref_fa_file, rebuild = False)
			self.chrlens = _np.array([len(x) for x in self.ref_fa_obj.records.values()]);
		except _FastaNotFoundError:
			print("Could not load reference genome!")

# export public functions
_fa = _FA(ref_fa_file = _os.environ["CAPY_REF_FA"] if "CAPY_REF_FA" in _os.environ else "/mnt/j/db/hg19/ref/hs37d5.fa");
genome_region = _fa._genome_region
chrpos2gpos = _fa._chrpos2gpos
gpos2chrpos = _fa._gpos2chrpos
get_chrlens = _fa._get_chrlens
set_reference = _fa._set_reference

#
# gnomAD functions
# 

 # TODO: update bin_stem to final when we generate with chrX
class _gnomad:
	def __init__(self, gnomad_dir = "/mnt/j/db/hg19/gnomad", bin_stem = "chr1-22", ref = None):
		self.gnomad_dir = gnomad_dir
		self.bin_stem = bin_stem
		self.mm_1bit = None
		self.ref = ref

		gnomad_index = gnomad_dir + "/1bit/" + bin_stem + ".index.parquet"
		if not _os.path.isfile(gnomad_index):
			print("Cannot find path to gnomAD index; gnomAD functionality disabled.", file = _sys.stderr)
			return

		self.obit_idx = _pd.read_parquet(gnomad_dir + "/1bit/" + bin_stem + ".index.parquet")

		self.__mmap()

	def __del__(self):
		if self.mm_1bit is not None:
			self.mm_1bit.close()

	def __mmap(self):
		if self.mm_1bit is not None:
			self.mm_1bit.close()

		f = open(self.gnomad_dir + "/1bit/" + self.bin_stem + ".bin", "rb");
		self.mm_1bit = _mmap.mmap(f.fileno(), 0, _mmap.MAP_SHARED, prot = _mmap.PROT_READ);
		f.close()

	def _query_1bit(self, ch, start, end = None):
		if end is None:
			start = start if not _np.isscalar(start) else [start]
			gpos_off = [(x,
			             *self.obit_idx.loc[
			               _np.digitize(x, self.obit_idx["g_start"]) - 1, ["offset", "g_start"]
			             ]
			            ) for x in chrpos2gpos(ch, start)
			           ]

			return [self.mm_1bit[offset + int((pos - g_start)/8)] &
			        (0x80 >> ((pos - g_start) % 8)) > 0 for pos, offset, g_start in gpos_off]
		else:
			pass
			#np.flatnonzero(np.unpackbits(np.frombuffer(self.mm_1bit[int(start/8):int(end/8)], dtype = np.uint8))) + 1

	def _query_1bit_raw(self, ch, start = None, end = None):
		if _np.size(_np.unique(ch)) > 1:
			raise ValueError("Cannot (yet) query multiple chromosomes simultaneously")

		if start is None:
			start = 1
		if end is None:
			end = get_chrlens(ref = self.ref)[ch - 1]

		offset = self.obit_idx.loc[_np.digitize(chrpos2gpos(ch, start), self.obit_idx["g_start"]) - 1, "offset"]

		return _np.frombuffer(self.mm_1bit[(offset + int(start/8)):(offset + int(_np.ceil(end/8)))], dtype = _np.uint8)

	# TODO: make this a decorator a la FASTA class?
	def _set_gnomad_ref_params(self, **kwargs):
		p = ["gnomad_dir", "bin_stem", "ref"]
		for arg, val in kwargs.items():
			if arg in p and val is not None:
				setattr(self, arg, val)

		self.__mmap()

_gnmd = _gnomad()
query_gnomad_1bit = _gnmd._query_1bit
query_gnomad_1bit_raw = _gnmd._query_1bit_raw
set_gnomad_ref_params = _gnmd._set_gnomad_ref_params
