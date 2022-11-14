import argparse
import copy
import multiprocessing
import pathlib
import pickle
import pysam
import sys
import time

import numpy as np

from source.tg_kmer import get_telomere_base_count, get_telomere_composition, read_kmer_tsv
from source.tg_plot import tel_len_violin_plot
from source.tg_tel  import get_double_anchored_tels, get_tels_below_canonical_thresh, parallel_anchored_tel_job
from source.tg_util import cluster_list, exists_and_is_nonzero, LEXICO_2_IND, makedir, parse_read, RC, repeated_matches_trimming

# hardcoded parameters
TEL_WINDOW_SIZE  = 100
# toss reads if too much of the non-telomere sequence couldn't be aligned anywhere
MAXIMUM_UNEXPLAINED_FRAC = 0.7
# toss reads if the median |p_vs_q| of the non-telomere sequence is above this value
MAX_NONTEL_MEDIAN_KMER_DENSITY = 0.25
#
MIN_DOUBLE_ANCHOR_LEN   = 1000
MIN_DOUBLE_ANCHOR_READS = 3
MIN_CANONICAL_FRAC      = 0.15

# for debugging purposes
DO_NOT_OVERWRITE = True

#
# ANCHORING_STRATEGY = 'largest'     - anchor tels onto largest non-tel alignment
# ANCHORING_STRATEGY = 'closest'     - anchor tels onto nearest non-tel alignment
#
# MATCH_TRIM_STRATEGY = 'largest'    - prioritize largest alignments when trimming overlaps
# MATCH_TRIM_STRATEGY = 'mapq'       - prioritize alignments with highest MAPQ when trimming overlaps
# MATCH_TRIM_STRATEGY = 'none'       - do not trim overlaps (use at your own risk for mm2 BAMs)
#

def main(raw_args=None):
	parser = argparse.ArgumentParser(description='Telogator v2.0', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
	parser.add_argument('-i',  type=str,   required=True,  metavar='input.bam',                      help="* Long reads aligned to t2t-telogator-ref.fa")
	parser.add_argument('-o',  type=str,   required=True,  metavar='output/',                        help="* Path to output directory")
	parser.add_argument('-r',  type=str,   required=False, metavar='hifi',        default='hifi',    help="Read type: hifi / clr / ont")
	parser.add_argument('-k',  type=str,   required=False, metavar='kmers.tsv',   default='',        help="Telomere kmers to search for")
	parser.add_argument('-t',  type=float, required=False, metavar='0.5',         default=0.5,       help="Telomere signal threshold (0-1)")
	#
	parser.add_argument('-fl', type=int,   required=False, metavar='5000',        default=5000,      help="Initial read filtering: Min read length")
	parser.add_argument('-ft', type=int,   required=False, metavar='200',         default=200,       help="Initial read filtering: Min tel bases")
	parser.add_argument('-ff', type=float, required=False, metavar='0.9',         default=0.9,       help="Anchored-tel filtering: Max frac of read that can be tel")
	parser.add_argument('-fp', type=float, required=False, metavar='0.25',        default=0.25,      help="Anchored-tel filtering: Max minor p/q fraction in tels")
	#
	parser.add_argument('-sa', type=str,   required=False, metavar='largest',     default='largest', help="Subtel/tel anchoring strategy")
	parser.add_argument('-sm', type=str,   required=False, metavar='mapq',        default='mapq',    help="Repeated matches trimming strategy")
	#
	parser.add_argument('-cd', type=int,   required=False, metavar='2000',        default=2000,      help="Maximum subtel distance to cluster anchored tels together")
	parser.add_argument('-cr', type=int,   required=False, metavar='3',           default=3,         help="Minimum number of reads per anchored-tel cluster")
	parser.add_argument('-ct', type=str,   required=False, metavar='p90',         default='p90',     help="Method for computing chr TL: mean / median / max / p90")
	#
	parser.add_argument('-vtm', type=int,  required=False, metavar='20000',       default=20000,     help="Violin plot (tel length) max value")
	parser.add_argument('-vtt', type=int,  required=False, metavar='5000',        default=5000,      help="Violin plot (tel length) tick size")
	parser.add_argument('-vrm', type=int,  required=False, metavar='50000',       default=50000,     help="Violin plot (read length) max value")
	parser.add_argument('-vrt', type=int,  required=False, metavar='10000',       default=10000,     help="Violin plot (read length) tick size")
	#
	parser.add_argument('--plot',          required=False, action='store_true',   default=False,     help="Create read plots")
	parser.add_argument('--plot-filt',     required=False, action='store_true',   default=False,     help="Create read plots (filtered reads)")
	parser.add_argument('--debug',         required=False, action='store_true',   default=False,     help="Print results for each read as its processed")
	#
	parser.add_argument('-p',   type=int,  required=False, metavar='4',           default=4,         help="Number of processes to use")

	args = parser.parse_args()

	#
	# initiate input / output files
	#
	INPUT_ALN  = args.i
	if exists_and_is_nonzero(INPUT_ALN) == False:
		print('Error reading -i input file')
		exit(1)
	if INPUT_ALN[-4:].lower() == '.sam':
		INPUT_TYPE = 'sam'
	elif INPUT_ALN[-4:].lower() == '.bam':
		INPUT_TYPE = 'bam'
	elif INPUT_ALN[-4:].lower() == '.cram':
		INPUT_TYPE = 'cram'
	elif INPUT_ALN[-2:].lower() == '.p':
		INPUT_TYPE = 'pickle'
	else:
		print('Error: unknown -i file type')
		exit(1)
	#
	OUT_DIR = args.o
	if OUT_DIR[-1] != '/':
		OUT_DIR += '/'
	makedir(OUT_DIR)
	OUT_PICKLE    = OUT_DIR + 'tel-data.p'
	OUT_PLOT_DIR  = OUT_DIR + 'read_plots/'
	OUT_CHR_TL    = OUT_DIR + 'tlens_by_chr.tsv'
	OUT_VIOLIN_TL = OUT_DIR + 'violin_tlens.png'
	OUT_VIOLIN_RL = OUT_DIR + 'violin_rlens.png'
	#
	PLOT_READS      = args.plot
	PLOT_FILT_READS = args.plot_filt
	PRINT_DEBUG     = args.debug
	#
	if PLOT_READS:
		makedir(OUT_PLOT_DIR)

	#
	# parse telomere kmer data
	#
	KMER_FILE = args.k
	if KMER_FILE != '':
		fn_suffix = KMER_FILE.split('/')[-1]
		print('using user-specified kmer list:', fn_suffix)
	else:
		print('using default telomere kmers')
		sim_path  = pathlib.Path(__file__).resolve().parent
		KMER_FILE = str(sim_path) + '/resources/kmers.tsv'
	(KMER_METADATA, KMER_ISSUBSTRING, CANONICAL_STRING) = read_kmer_tsv(KMER_FILE)
	[KMER_LIST, KMER_COLORS, KMER_LETTER, KMER_FLAGS] = KMER_METADATA
	KMER_LIST_REV = [RC(n) for n in KMER_LIST]
	CANONICAL_STRING_REV = RC(CANONICAL_STRING)

	#
	# various parameters
	#
	READ_TYPE         = args.r
	MINIMUM_READ_LEN  = args.fl
	MINIMUM_TEL_BASES = args.ft
	MAXIMUM_TEL_FRAC  = args.ff
	MAXIMUM_MINOR_PQ  = args.fp
	P_VS_Q_AMP_THRESH = args.t
	#
	ANCHORING_STRATEGY  = args.sa
	MATCH_TRIM_STRATEGY = args.sm
	#
	ANCHOR_CLUSTER_DIST = args.cd
	MIN_READS_PER_CLUST = args.cr
	CHR_TL_METHOD       = args.ct
	#
	(VIOLIN_TLEN_MAX, VIOLIN_TLEN_TICK) = (args.vtm, args.vtt)
	(VIOLIN_RLEN_MAX, VIOLIN_RLEN_TICK) = (args.vrm, args.vrt)
	#
	NUM_PROCESSES = args.p

	reads_skipped = {'min_readlen':0,
	                 'trim_filter':0,
	                 'min_telbases':0}

	ANCHORED_TEL_BY_CHR    = {}
	NONTEL_REFSPANS_BY_CHR = {}

	#
	# GET ALL READ DATA FROM INPUT ALIGNMENT
	#
	ALIGNMENTS_BY_RNAME = {}
	#
	if INPUT_TYPE == 'pickle':
		print('reading pickle input...')
		f = open(INPUT_ALN, 'rb')
		my_pickle = pickle.load(f)
		f.close()
		ANCHORED_TEL_BY_CHR = my_pickle['anchored-tels']
		NONTEL_REFSPANS_BY_CHR   = my_pickle['non-tel-ref-spans']
	else:
		if INPUT_TYPE == 'sam':
			samfile = pysam.AlignmentFile(INPUT_ALN, "r")
		elif INPUT_TYPE == 'bam':
			samfile = pysam.AlignmentFile(INPUT_ALN, "rb")
		elif INPUT_TYPE == 'cram':
			samfile = pysam.AlignmentFile(INPUT_ALN, "rc")
		#
		refseqs = samfile.references
		sys.stdout.write('getting all read data from input alignment...')
		sys.stdout.flush()
		tt = time.time()
		for aln in samfile.fetch(until_eof=True):
			sam_line    = str(aln).split('\t')
			sam_line[2] = refseqs[int(sam_line[2])]	# pysam is dumb and prints ref indices instead of contig name
			#
			[rnm, ref_key, pos, read_pos_1, read_pos_2, ref, pos1, pos2, orientation, mapq, rdat] = parse_read(sam_line)
			if rnm not in ALIGNMENTS_BY_RNAME:
				ALIGNMENTS_BY_RNAME[rnm] = []
			ALIGNMENTS_BY_RNAME[rnm].append([read_pos_1, read_pos_2, ref, pos1, pos2, orientation, mapq, rdat])
		sys.stdout.write(' (' + str(int(time.time() - tt)) + ' sec)\n')
		sys.stdout.flush()

		#
		# INITIAL FILTERING
		#
		FILTERED_READS = []
		#
		sys.stdout.write('initial read filtering...')
		sys.stdout.flush()
		num_starting_reads = len(ALIGNMENTS_BY_RNAME.keys())
		tt = time.time()
		for readname in ALIGNMENTS_BY_RNAME.keys():
			abns_k = repeated_matches_trimming(sorted(ALIGNMENTS_BY_RNAME[readname]), strategy=MATCH_TRIM_STRATEGY, print_debug=PRINT_DEBUG)
			# did we lose all of our alignments during trimming?
			if len(abns_k) == 0:
				reads_skipped['trim_filter'] += 1
				continue
			# make sure string used for kmer matching is same orientation as the alignments
			which_i = 0
			for i in range(len(abns_k)):
				if abns_k[i][2][:3] != 'tel':
					which_i = i
					break
			# assuming softclipping was used. i.e. all alignments should have same sequence... (don't pick tel though)
			if abns_k[which_i][5] == 'FWD':
				rdat = abns_k[which_i][7]
			elif abns_k[which_i][5] == 'REV':
				rdat = RC(abns_k[which_i][7])
			# read len filter
			if len(rdat) < MINIMUM_READ_LEN:
				reads_skipped['min_readlen'] += 1
				continue
			# check if we're unmapped
			refs_we_aln_to = [aln[2] for aln in abns_k]
			refs_we_aln_to = sorted(list(set(refs_we_aln_to)))
			if refs_we_aln_to == ['*']:
				reads_skipped['unmapped'] += 1
				continue
			# minimum tel content
			tel_bc = get_telomere_base_count(rdat, [CANONICAL_STRING, CANONICAL_STRING_REV], mode=READ_TYPE)
			if tel_bc < MINIMUM_TEL_BASES:
				reads_skipped['min_telbases'] += 1
				if INPUT_TYPE != 'pickle':	# use non-tel reads for downstream filtering of double-anchored tels
					for aln in abns_k:
						if aln[2][:3] != 'tel':
							if aln[2] not in NONTEL_REFSPANS_BY_CHR:
								NONTEL_REFSPANS_BY_CHR[aln[2]] = []
							NONTEL_REFSPANS_BY_CHR[aln[2]].append(tuple(sorted(aln[3:5])))
				continue
			# we passed all filters?
			FILTERED_READS.append([readname, rdat, copy.deepcopy(abns_k)])
		#
		sys.stdout.write(' (' + str(int(time.time() - tt)) + ' sec)\n')
		num_ending_reads = len(FILTERED_READS)
		sys.stdout.write(' - ' + str(num_starting_reads) + ' --> ' + str(num_ending_reads) + ' reads\n')
		sys.stdout.flush()

		#
		#
		# TELOGATOR 1.0 - IDENTIFY ANCHORED TELOMERES
		#
		#
		if NUM_PROCESSES <= 1:
			all_indices      = list(range(len(FILTERED_READS)))
			reads_to_process = FILTERED_READS
		else:
			all_indices      = [[] for n in range(NUM_PROCESSES)]
			reads_to_process = [[] for n in range(NUM_PROCESSES)]
			k = 0
			for i in range(len(FILTERED_READS)):
				all_indices[k%NUM_PROCESSES].append(i)
				reads_to_process[k%NUM_PROCESSES].append(copy.deepcopy(FILTERED_READS[i]))
				k += 1
		#
		# execute parallel jobs
		#
		sys.stdout.write('identify anchorable telomeres...')
		sys.stdout.flush()
		num_starting_reads = len(FILTERED_READS)
		num_ending_reads   = 0
		tt = time.time()
		#
		par_params      = [KMER_LIST, KMER_LIST_REV, TEL_WINDOW_SIZE, P_VS_Q_AMP_THRESH, ANCHORING_STRATEGY, PLOT_READS, INPUT_TYPE, OUT_PLOT_DIR, PRINT_DEBUG, PLOT_FILT_READS]
		par_params_filt = [MAXIMUM_TEL_FRAC, MAXIMUM_MINOR_PQ, MAXIMUM_UNEXPLAINED_FRAC, MAX_NONTEL_MEDIAN_KMER_DENSITY]
		if NUM_PROCESSES <= 1:
			par_results = {}
			parallel_anchored_tel_job(reads_to_process, all_indices, par_params, par_params_filt, par_results)
		else:
			manager     = multiprocessing.Manager()
			par_results = manager.dict()
			processes   = []
			for i in range(NUM_PROCESSES):
				p = multiprocessing.Process(target=parallel_anchored_tel_job, args=(reads_to_process[i], all_indices[i], par_params, par_params_filt, par_results))
				processes.append(p)
			for i in range(NUM_PROCESSES):
				processes[i].start()
			for i in range(NUM_PROCESSES):
				processes[i].join()
		#
		# collect parallel job results
		#
		for i in range(k):
			anchored_tel_dat = par_results[(i,0)]
			my_filt_string   = par_results[(i,1)]
			nontel_spans     = par_results[(i,2)]
			#
			if len(my_filt_string):
				if my_filt_string not in reads_skipped:
					reads_skipped[my_filt_string] = 0
				reads_skipped[my_filt_string] += 1
			else:
				adjacent_chr = anchored_tel_dat[0]
				if adjacent_chr not in ANCHORED_TEL_BY_CHR:
					ANCHORED_TEL_BY_CHR[adjacent_chr] = []
				ANCHORED_TEL_BY_CHR[adjacent_chr].append(anchored_tel_dat[1:])
				num_ending_reads += 1
			#
			for k in nontel_spans.keys():
				if k not in NONTEL_REFSPANS_BY_CHR:
					NONTEL_REFSPANS_BY_CHR[k] = []
				NONTEL_REFSPANS_BY_CHR[k].extend(nontel_spans[k])
		#
		sys.stdout.write(' (' + str(int(time.time() - tt)) + ' sec)\n')
		sys.stdout.write(' - ' + str(num_starting_reads) + ' --> ' + str(num_ending_reads) + ' reads\n')
		sys.stdout.flush()

		#
		# sort non-tel ref spans
		#
		for k in NONTEL_REFSPANS_BY_CHR.keys():
			NONTEL_REFSPANS_BY_CHR[k] = sorted(NONTEL_REFSPANS_BY_CHR[k])
		#
		# print stats on reads that were filtered
		#
		print()
		for k in sorted(reads_skipped.keys()):
			if reads_skipped[k] > 0:
				print(k, reads_skipped[k])
		print()
		#
		# save output pickle
		#
		f = open(OUT_PICKLE, 'wb')
		pickle.dump({'anchored-tels':ANCHORED_TEL_BY_CHR, 'non-tel-ref-spans':NONTEL_REFSPANS_BY_CHR}, f)
		f.close()

	#
	# remove telomeres that we thought were single-anchored but other reads tell us it is actually double-anchored
	#
	sys.stdout.write('removing anchored tels that double-anchored according to other reads...')
	sys.stdout.flush()
	num_starting_reads = sum([len(ANCHORED_TEL_BY_CHR[k]) for k in ANCHORED_TEL_BY_CHR.keys()])
	tt = time.time()
	#
	gdat_params = [MIN_DOUBLE_ANCHOR_LEN, MIN_DOUBLE_ANCHOR_READS, PRINT_DEBUG]
	del_keys    = get_double_anchored_tels(ANCHORED_TEL_BY_CHR, NONTEL_REFSPANS_BY_CHR, gdat_params)
	#
	num_ending_reads = num_starting_reads - len(del_keys)
	del_keys2 = []
	for (k,di) in del_keys:
		del ANCHORED_TEL_BY_CHR[k][di]
		if len(ANCHORED_TEL_BY_CHR[k]) == 0:
			del_keys2.append(k)
	for k in del_keys2:
		del ANCHORED_TEL_BY_CHR[k]
	sys.stdout.write(' (' + str(int(time.time() - tt)) + ' sec)\n')
	sys.stdout.write(' - ' + str(num_starting_reads) + ' --> ' + str(num_ending_reads) + ' reads\n')
	sys.stdout.flush()

	#
	# remove reads where telomere region is insufficiently canonical
	#
	sys.stdout.write('removing reads with tel regions that are insufficiently canonical...')
	sys.stdout.flush()
	num_starting_reads = sum([len(ANCHORED_TEL_BY_CHR[k]) for k in ANCHORED_TEL_BY_CHR.keys()])
	#
	gtbct_params = [MIN_CANONICAL_FRAC, CANONICAL_STRING, CANONICAL_STRING_REV, READ_TYPE]
	del_keys     = get_tels_below_canonical_thresh(ANCHORED_TEL_BY_CHR, gtbct_params)
	#
	num_ending_reads = num_starting_reads - len(del_keys)
	del_keys2 = []
	for (k,di) in del_keys:
		del ANCHORED_TEL_BY_CHR[k][di]
		if len(ANCHORED_TEL_BY_CHR[k]) == 0:
			del_keys2.append(k)
	for k in del_keys2:
		del ANCHORED_TEL_BY_CHR[k]
	sys.stdout.write(' (' + str(int(time.time() - tt)) + ' sec)\n')
	sys.stdout.write(' - ' + str(num_starting_reads) + ' --> ' + str(num_ending_reads) + ' reads\n')
	sys.stdout.flush()

	#
	# cluster subtel/tel boundaries by position and compute TL statistics
	#
	sorted_ref_keys = []
	for k in sorted(ANCHORED_TEL_BY_CHR.keys()):
		if k[:3] == 'chr':
			my_chr = k.replace('_','')
			my_ind = 1
		elif k[:3] == 'alt':
			my_chr = ''.join(k.split('_')[:-1])
			my_chr = my_chr.replace('alt', 'chr')
			my_ind = 2
		else:
			print('skipping telomeres anchored to unknown ref:', k)
			continue
		sorted_ref_keys.append((my_ind, LEXICO_2_IND[my_chr[:-1]], my_chr, k))
	sorted_ref_keys = sorted(sorted_ref_keys)
	#
	print('clustering anchored tels by position & filtering by read count...')
	CHR_TEL_DAT = []
	TEL_LEN_OUT = {}	# indexed by chr
	READLEN_OUT = {}	# indexed by chr
	clust_num   = 0
	tel_composition_data = []
	for ki in range(len(sorted_ref_keys)):
		(my_chr, k) = (sorted_ref_keys[ki][2], sorted_ref_keys[ki][3])
		sort_list = sorted([(ANCHORED_TEL_BY_CHR[k][n][1], n) for n in range(len(ANCHORED_TEL_BY_CHR[k]))])
		clusters  = cluster_list(sort_list, ANCHOR_CLUSTER_DIST, which_val=0)
		num_reads = sum([len(n) for n in clusters])
		sys.stdout.write(' - ' + my_chr + ': ' + str(len(clusters)) + ' anchor / ' + str(num_reads) + ' read')
		num_clust = 0
		num_reads = 0
		for cl in clusters:
			if len(cl) < MIN_READS_PER_CLUST:
				continue
			num_clust += 1
			num_reads += len(cl)
			pos_list = [n[0] for n in cl]
			ind_list = [n[1] for n in cl]
			my_pos   = int(np.median(pos_list))
			clust_num += 1
			#
			sorted_by_tlen = sorted([(ANCHORED_TEL_BY_CHR[k][i][3], len(ANCHORED_TEL_BY_CHR[k][i][6]), ANCHORED_TEL_BY_CHR[k][i][0]) for i in ind_list])
			my_tlens       = [n_sbtl[0] for n_sbtl in sorted_by_tlen]
			my_rlens       = [n_sbtl[1] for n_sbtl in sorted_by_tlen]
			my_rnames      = [n_sbtl[2] for n_sbtl in sorted_by_tlen]
			#
			if CHR_TL_METHOD == 'mean':
				consensus_tl = np.mean(my_tlens)
			elif CHR_TL_METHOD == 'median':
				consensus_tl = np.median(my_tlens)
			elif CHR_TL_METHOD == 'max':
				consensus_tl = np.max(my_tlens)
			elif CHR_TL_METHOD[0] == 'p':
				my_percentile = int(CHR_TL_METHOD[1:])
				consensus_tl = np.percentile(my_tlens, my_percentile)
			else:
				print('Error: unknown chr TL summary method')
				exit(1)
			CHR_TEL_DAT.append([my_chr,
			                    str(my_pos),
			                    str(int(consensus_tl)),
			                    ','.join([str(n) for n in my_tlens]),
			                    ','.join([str(n) for n in my_rlens])])
			#
			if my_chr not in TEL_LEN_OUT:
				TEL_LEN_OUT[my_chr] = []
				READLEN_OUT[my_chr] = []
			TEL_LEN_OUT[my_chr].extend(my_tlens)
			READLEN_OUT[my_chr].extend(my_rlens)
			#
			clust_tcdat = []
			gtc_params  = [my_chr, clust_num, KMER_LIST, KMER_LIST_REV, KMER_ISSUBSTRING]
			for i in ind_list:
				clust_tcdat.append(get_telomere_composition(ANCHORED_TEL_BY_CHR[k][i], gtc_params))
			clust_tcout = [my_chr, clust_num, ind_list, clust_tcdat]
			tel_composition_data.append(copy.deepcopy(clust_tcout))
			#
		if num_clust != len(clusters):
			sys.stdout.write(' --> ' + str(num_clust) + ' anchor / ' + str(num_reads) + ' read\n')
		else:
			sys.stdout.write('\n')
		sys.stdout.flush()

	#
	# write out chromosome TLs
	#
	print('writing TL results to "' + OUT_CHR_TL.split('/')[-1] + '"...')
	f = open(OUT_CHR_TL, 'w')
	f.write('#chr' + '\t' +
		    'position' + '\t' +
		    'TL_' + CHR_TL_METHOD + '\t' +
		    'read_TLs' + '\t' +
		    'read_lengths' + '\n')
	for i in range(len(CHR_TEL_DAT)):
		f.write('\t'.join(CHR_TEL_DAT[i]) + '\n')
	f.close()

	#
	# violin plots
	#
	print('creating violin plots...')
	tlen_violin_params = {'p_ymax':VIOLIN_TLEN_MAX,
	                      'q_ymax':VIOLIN_TLEN_MAX,
	                      'y_step':VIOLIN_TLEN_TICK}
	rlen_violin_params = {'p_color':'gray',
	                      'q_color':'gray',
	                      'y_label':'<-- q     read length     p -->',
	                      'p_ymax':VIOLIN_RLEN_MAX,
	                      'q_ymax':VIOLIN_RLEN_MAX,
	                      'y_step':VIOLIN_RLEN_TICK}
	tel_len_violin_plot(TEL_LEN_OUT, OUT_VIOLIN_TL, custom_plot_params=tlen_violin_params)
	tel_len_violin_plot(READLEN_OUT, OUT_VIOLIN_RL, custom_plot_params=rlen_violin_params)

if __name__ == '__main__':
	main()
