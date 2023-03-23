import argparse
import gzip

from source.tg_reader import TG_Reader
from source.tg_util   import RC

def main(raw_args=None):
	parser = argparse.ArgumentParser(description='get_kmer_hits.py', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
	parser.add_argument('-i', type=str, required=True,  metavar='input.fa',                    help="* Input reads (.fa / .fa.gz / .fq / .fq.gz)")
	parser.add_argument('-o', type=str, required=True,  metavar='output.fa',                   help="* Output reads (.fa / .fa.gz / .fq / .fq.gz)")
	parser.add_argument('-k', type=str, required=False, metavar='TTAGGG',    default='TTAGGG', help="Kmer to search")
	parser.add_argument('-m', type=int, required=False, metavar='20',        default=20,       help="Minimum hit count to output read")
	args = parser.parse_args()
	#
	IN_READS  = args.i
	OUT_READS = args.o
	KMER      = args.k
	KMER_RC   = RC(KMER)
	MIN_HITS  = args.m
	#
	OUTPUT_IS_FASTQ = False
	if OUT_READS[-3:].lower() == '.fq' or OUT_READS[-6:].lower() == '.fq.gz' or OUT_READS[-6:] == '.fastq' or OUT_READS[-9:] == '.fastq.gz':
		OUTPUT_IS_FASTQ = True
	INPUT_IS_FASTQ = False
	if IN_READS[-3:].lower() == '.fq' or IN_READS[-6:].lower() == '.fq.gz' or IN_READS[-6:] == '.fastq' or IN_READS[-9:] == '.fastq.gz':
		INPUT_IS_FASTQ = True
	if INPUT_IS_FASTQ == False and OUTPUT_IS_FASTQ == True:
		print('Error: input is fasta and output is fastq.')
		exit(1)
	#
	if OUT_READS[-3:].lower() == '.gz':
		f_out = gzip.open(OUT_READS, 'wt')
	else:
		f_out = open(OUT_READS, 'w')
	my_reader = TG_Reader(IN_READS)
	#
	hit_count_dict = {}
	while True:
		(my_name, my_rdat, my_qdat) = my_reader.get_next_read()
		if not my_name:
			break
		count_fwd = my_rdat.count(KMER)
		count_rev = my_rdat.count(KMER_RC)
		if count_fwd not in hit_count_dict:
			hit_count_dict[count_fwd] = 0
		hit_count_dict[count_fwd] += 1
		if count_rev not in hit_count_dict:
			hit_count_dict[count_rev] = 0
		hit_count_dict[count_rev] += 1
		if count_fwd > MIN_HITS or count_rev > MIN_HITS:
			if OUTPUT_IS_FASTQ:
				f_out.write('@' + my_name + '\n' + my_rdat + '\n' + '+' + '\n' + my_qdat + '\n')
			else:
				f_out.write('>' + my_name + '\n' + my_rdat + '\n')
	#
	f_out.close()
	my_reader.close()
	#
	for k in sorted(hit_count_dict.keys()):
		print(k, hit_count_dict[k])

if __name__ == '__main__':
	main()
