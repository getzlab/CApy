from pyfaidx import Fasta as _Fasta
import numpy as _np
import sys as _sys
import mmap as _mmap
import pandas as _pd

# wrap in a class because we only want to remember which reference we're using,
# so we don't have to instantiate it every time genome_region is called
class _FA:
	def __init__(self, ref_fa_file = "/mnt/j/db/hg19/ref/hs37d5.fa"):
		self.__set_reference(ref_fa_file)

	def __auto_reset_reference(f):
		def inner(self, *args, ref = None, **kwargs):
			if ref is None:
				ref = self.ref_fa_file

			if ref != self.ref_fa_file:
				self.__set_reference(ref)
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
	def _get_chrlens(self, ref = None):
		return self.chrlens

	def __set_reference(self, ref):
		self.ref_fa_file = ref
		self.ref_fa_obj = _Fasta(self.ref_fa_file)
		self.chrlens = _np.array([len(x) for x in self.ref_fa_obj.records.values()]);

# export public functions
_fa = _FA();
genome_region = _fa._genome_region
chrpos2gpos = _fa._chrpos2gpos
get_chrlens = _fa._get_chrlens

class _gnomad:
	def __init__(self, gnomad_dir = "/mnt/j/db/hg19/gnomad", ref = None):
		self.gnomad_dir = gnomad_dir
		self.obit_idx = _pd.read_parquet(gnomad_dir + "/1bit/chr1-22.index.parquet") # TODO: update to final when we generate with chrX
		self.ref = ref

		f = open(self.gnomad_dir + "/1bit/chr1-22.bin", "rb"); # TODO: update to final when we generate with chrX
		self.mm_1bit = _mmap.mmap(f.fileno(), 0, _mmap.MAP_SHARED, prot = _mmap.PROT_READ);
		f.close()

	def __del__(self):
		self.mm_1bit.close()

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

#	def _query_1bit_raw(self, ch, start, end):
#		[gst_off, gen_off] = [(x, np.digitize(x, self.obit_idx["g_start"]) - 1) for x in chrpos2gpos(ch, [start, end])]
#
#		return np.frombuffer(self.mm_1bit[int(start/8):int(end/8)], dtype = np.uint8)

_gnmd = _gnomad()
query_gnomad_1bit = _gnmd._query_1bit
