from pyfaidx import Fasta as _Fasta
import numpy as _np
import sys as _sys

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

	def __set_reference(self, ref):
		self.ref_fa_file = ref
		self.ref_fa_obj = _Fasta(self.ref_fa_file)
		self.chrlens = _np.array([len(x) for x in self.ref_fa_obj.records.values()]);

# export public functions
_fa = _FA();
genome_region = _fa._genome_region
chrpos2gpos = _fa._chrpos2gpos
