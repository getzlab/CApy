from pyfaidx import Fasta
import sys

# wrap in a class because we only want to remember which reference we're using,
# so we don't have to instantiate it every time genome_region is called
class FA:
	def __init__(self, ref_fa_file = "/mnt/j/db/hg19/ref/hs37d5.fa"):
		self.ref_fa_file = ref_fa_file
		self.ref_fa_obj = Fasta(self.ref_fa_file)

	def _genome_region(self, chr, start, end, ref = None):
		if ref is None:
			ref = self.ref_fa_file

		if ref != self.ref_fa_file:
			self.__set_reference(ref)
			print("NOTE: switched reference to {}".format(ref), file = sys.stderr)

		# TODO: handle implicit chr conversion better
		return self.ref_fa_obj[str(chr)][start:end].seq

	def __set_reference(self, ref):
		self.ref_fa_file = ref
		self.ref_fa_obj = Fasta(self.ref_fa_file)

# export public functions
fa = FA();
genome_region = fa._genome_region
