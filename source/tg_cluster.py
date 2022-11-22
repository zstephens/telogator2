import os
import multiprocessing
import random

import numpy as np
import matplotlib.pyplot as mpl

from collections import Counter

from Bio import pairwise2

from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist, squareform

from source.tg_reader import TG_Reader
from source.tg_util import exists_and_is_nonzero

# we're going to pretend each kmer color is an amino acid, so alignment tools behave themselves
AMINO = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
AMINO_2_IND = {AMINO[i]:i for i in range(len(AMINO))}
UNKNOWN_LETTER = AMINO[0]

# how many blocks of gaps to consider for plotting adj
MAX_GAP_BLOCK = 3

# how many shuffles to perform for alignment score background model
RAND_SHUFFLE_COUNT = 3

# when comparing sequences of different lengths, choose min(len(seq1), len(seq2)), but only go this low
MIN_VIABLE_SEQ_LEN = 1000
# the log-based distance function can go to infinity so lets set an upper bound on it
MAX_SEQ_DIST = 10.0
# similarly, lets choose a small number as the minimum to prevent numerical weirdness from giving us negative values
MIN_SEQ_DIST = 0.0001
# to prevent pesky div-by-zeros in edge cases
MIN_MSD      = 3.0

MATCH_NORMAL  = 5
XMATCH_NORMAL = -4
GAP_OPEN = -4
GAP_EXT  = -4
#
MATCH_CANON  = 0
XMATCH_CANON = -4
#
MATCH_UNKNOWN  = 2
XMATCH_UNKNOWN = -4

DEFAULT_TREECUT = 0.250

# density parameters for identifing subtel / tvr boundaries
UNKNOWN_WIN_SIZE = 100
UNKNOWN_END_DENS = 0.120
# density parameters for discerning canonical regions from sequencing artifacts
CANON_WIN_SIZE = 100
CANON_END_DENS = 0.700
# parameters for determining tvr / canonical boundaries: (denoise_region_size, cum_thresh, min_hits)
TVR_CANON_FILT_PARAMS_STRICT  = (10, 0.05, 100)
TVR_CANON_FILT_PARAMS_LENIENT = ( 5, 0.20,  50)
#
MAX_TVR_LEN       = 8000	# ignore variant repeats past this point when finding TVR boundary
MAX_TVR_LEN_SHORT = 3000	# when examining TVRs with very few variant repeats

def write_amino_scoring_matrix(fn):
	f = open(fn, 'w')
	f.write('   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *' + '\n')
	f.write('A  1  0  0  0  5  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -4' + '\n')
	f.write('R  0  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('N  0 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('D  0 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('C  5 -4 -4 -4  0 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('Q  0 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('E  0 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('G  0 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('H  0 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('I  0 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('L  0 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('K  0 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('M  0 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('F  0 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('P  0 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('S  0 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('T  0 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('W  0 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('Y  0 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('V  0 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -1 -4' + '\n')
	f.write('B  0 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -1 -4' + '\n')
	f.write('J  0 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -1 -4' + '\n')
	f.write('Z  0 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -1 -4' + '\n')
	f.write('X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4' + '\n')
	f.write('* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1' + '\n')
	f.close()

def write_nucl_scoring_matrix(fn):
	f = open(fn, 'w')
	f.write('   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *' + '\n')
	f.write('A  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  0 -4' + '\n')
	f.write('R -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('N -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('D -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('C -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('Q -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('E -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('G -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('H -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('I -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('L -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('K -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('M -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('F -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('P -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('S -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('T -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('W -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('Y -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('V -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -1 -4' + '\n')
	f.write('B -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -1 -4' + '\n')
	f.write('J -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -1 -4' + '\n')
	f.write('Z -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -1 -4' + '\n')
	f.write('X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4' + '\n')
	f.write('* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1' + '\n')
	f.close()

def get_muscle_msa(input_sequences, muscle_exe, working_dir='', char_score_adj={}, max_gap_frac=0.60, mode='amino'):
	# write sequences to a temp fasta
	temp_fasta = working_dir + 'clust_sequences.fa'
	f = open(temp_fasta, 'w')
	for i in range(len(input_sequences)):
		f.write('>seq'+str(i+1).zfill(5) + '\n')
		f.write(input_sequences[i] + '\n')
	f.close()
	# run muscle
	aln_fasta   = working_dir + 'clust_aln.fa'
	muscle_log  = working_dir + 'muscle.log'
	matrix      = working_dir + 'scoring_matrix.txt'
	if mode == 'amino':
		write_amino_scoring_matrix(matrix)
		score_param = '-seqtype protein -gapopen -12.0 -gapextend -4.0 -center 0.0 -matrix ' + matrix
	elif mode == 'nucl':
		write_nucl_scoring_matrix(matrix)
		score_param = '-seqtype dna -gapopen -12.0 -gapextend -4.0 -center 0.0 -matrix ' + matrix
	else:
		print('Error: get_muscle_msa mode must be amino or nucl')
		exit(1)
	score_param += ' -maxiters 2'	# use if muscle is crashing on "refining bipartite" steps
	cmd = muscle_exe + ' -in ' + temp_fasta + ' -out ' + aln_fasta + ' ' + score_param + ' > ' + muscle_log + ' 2>&1'
	os.system(cmd)
	# get results
	out_seq   = []
	my_reader = TG_Reader(aln_fasta, verbose=False)
	while True:
		read_dat = my_reader.get_next_read()
		if not read_dat[0]:
			break
		my_readnum = int(read_dat[0][3:])
		out_seq.append([my_readnum, read_dat[1]])
	my_reader.close()
	out_seq = [n[1] for n in sorted(out_seq)]
	####for n in out_seq:
	####	print(n)
	####exit(1)
	# cleanup
	os.system('rm ' + temp_fasta)
	os.system('rm ' + aln_fasta)
	os.system('rm ' + muscle_log)
	os.system('rm ' + matrix)
	# get consensus
	consensus_seq = []
	for i in range(len(out_seq[0])):	# we're gonna hope that muscle worked as expected and all seq are same len
		char_count = {}
		gap_count  = 0
		for j in range(len(out_seq)):
			if out_seq[j][i] == '-':
				gap_count += 1
			else:
				if out_seq[j][i] not in char_count:
					char_count[out_seq[j][i]] = 0
				char_count[out_seq[j][i]] += 1
		#
		if float(gap_count)/len(out_seq) > max_gap_frac:
			continue
		#
		candidates = [(char_count[k],k) for k in char_count.keys() if char_count[k] == max(char_count.values())]
		if len(candidates) == 1:
			consensus_seq.append(candidates[0][1])
		else:	# tie-breaking logic
			adj_scores = []
			for candidate in candidates:
				if candidate[1] in char_score_adj:
					adj_scores.append((char_score_adj[candidate[1]], candidate[1]))
				else:
					adj_scores.append((0, candidate[1]))
			adj_scores = sorted(adj_scores, reverse=True)
			#print(candidates, '-->', adj_scores)
			consensus_seq.append(adj_scores[0][1])
	consensus_seq = ''.join(consensus_seq)
	#
	return [out_seq, consensus_seq]

#
#
#
def get_adj_from_gaps(s):
	in_gap     = False
	gap_count  = 0
	gap_blocks = []
	for i in range(len(s)):
		if s[i] == '-':
			if in_gap == False:
				gap_start = i
			in_gap = True
		else:
			if in_gap:
				gap_blocks.append((gap_start,i))
				gap_count += 1
			in_gap = False
		if gap_count >= MAX_GAP_BLOCK:
			break
	out_adj = 0
	if len(gap_blocks):
		gap_so_far  = [gap_blocks[0][1]-gap_blocks[0][0]]
		cost_so_far = [gap_blocks[0][0]]
		for i in range(1,len(gap_blocks)):
			gap_so_far.append(gap_so_far[-1] + gap_blocks[i][1]-gap_blocks[i][0])
			cost_so_far.append(cost_so_far[-1] + gap_blocks[i][0]-gap_blocks[i-1][1])
		best_gap = sorted([(gap_so_far[i]-cost_so_far[i],i) for i in range(len(gap_blocks))], reverse=True)[0]
		if best_gap[0] > 0:
			out_adj = gap_blocks[best_gap[1]][1]
	return out_adj

def shuffle_seq(s):
	return ''.join(random.sample(s,len(s)))

#####
##### linear regression parameters for estimating rand alignment score (old colors)
#####
####RAND_REG_BIAS    =  82.745688430256450
####RAND_REG_WEIGHTS = [ 0.318626543939361,
####                    -0.109363196930245,
####                     0.049521602321889,
####                    -0.546815984651203,
####                     0.058014516390758,
####                    -1.595058825106752,
####                     6.719760869407775]
#
# linear regression parameters for estimating rand alignment score (new colors)
#
RAND_REG_BIAS    = 156.660399273161600
RAND_REG_WEIGHTS = [ 0.360900872720303,
                    -0.115378786823871,
                     0.043744885084501,
                    -0.592517362874134,
                     0.221927915784329,
                    -2.547859906771088,
                     6.772960422330248]
#
def estimate_rand_score(seq1, seq2, aln_score, iden_score):
	c1 = Counter(seq1)
	c2 = Counter(seq2)
	for letter in ['A','C']:
		if letter not in c1:
			c1[letter] = 0
		if letter not in c2:
			c2[letter] = 0
	common_letters = c1 & c2
	all_common     = sum(common_letters.values())
	a_common       = min(c1['A'], c2['A'])
	c_common       = min(c1['C'], c2['C'])
	len_min        = min(len(seq1), len(seq2))
	len_max        = max(len(seq1), len(seq2))
	feature_vector = [len_min, len_max, aln_score, iden_score, all_common, a_common, c_common]
	rand_score = RAND_REG_BIAS
	for i in range(len(feature_vector)):
		rand_score += RAND_REG_WEIGHTS[i]*feature_vector[i]
	return int(rand_score)

#
#
#
def parallel_alignment_job(our_indices, sequences, gap_bool, pq, out_dict, scoring_matrix=None, train_out=None, estimate_rand=False):
	for (i,j) in our_indices:
		#
		min_len = min(len(sequences[i]), len(sequences[j]))
		min_len = max(min_len, MIN_VIABLE_SEQ_LEN)
		#
		if pq == 'p':
			seq_i = sequences[i][-min_len:]
			seq_j = sequences[j][-min_len:]
		elif pq == 'q':
			seq_i = sequences[i][:min_len]
			seq_j = sequences[j][:min_len]
		#
		# aln score
		#
		if scoring_matrix == None:
			aln = pairwise2.align.globalms(seq_i, seq_j, MATCH_NORMAL, XMATCH_NORMAL, GAP_OPEN, GAP_EXT, penalize_end_gaps=gap_bool)
		else:
			aln = pairwise2.align.globalds(seq_i, seq_j, scoring_matrix, GAP_OPEN, GAP_EXT, penalize_end_gaps=gap_bool)
		aln_score = int(aln[0].score)
		#
		# iden score
		#
		if scoring_matrix == None:
			iden_score = MATCH_NORMAL*max(len(seq_i), len(seq_j))
		else:
			is1 = sum([scoring_matrix[(n,n)] for n in seq_i])
			is2 = sum([scoring_matrix[(n,n)] for n in seq_j])
			iden_score = max(is1, is2)
		#
		# rand score
		#
		if estimate_rand:
			rand_score = estimate_rand_score(seq_i, seq_j, aln_score, iden_score)
		else:
			rand_scores = []
			for k in range(RAND_SHUFFLE_COUNT):
				if scoring_matrix == None:
					rand_aln = pairwise2.align.globalms(shuffle_seq(seq_i), shuffle_seq(seq_j), MATCH_NORMAL, XMATCH_NORMAL, GAP_OPEN, GAP_EXT, penalize_end_gaps=gap_bool)
				else:
					rand_aln = pairwise2.align.globalds(shuffle_seq(seq_i), shuffle_seq(seq_j), scoring_matrix, GAP_OPEN, GAP_EXT, penalize_end_gaps=gap_bool)
				rand_scores.append(rand_aln[0].score)
			rand_score = int(np.mean(rand_scores))
		#
		# distance calculation
		#
		if rand_score >= aln_score:
			my_dist = MAX_SEQ_DIST
		else:
			my_dist = min(-np.log((aln_score - rand_score)/(iden_score - rand_score)), MAX_SEQ_DIST)
		my_dist = max(my_dist, MIN_SEQ_DIST)
		#
		# output
		#
		print((i,j), aln_score, rand_score, iden_score, '{0:0.3f}'.format(my_dist))
		if train_out != None:
			f = open(train_out, 'a')
			f.write(seq_i + '\t' + seq_j + '\t' + str(aln_score) + '\t' + str(iden_score) + '\t' + str(rand_score) + '\n')
			f.close()
		out_dict[(i,j)] = my_dist

#
#
#
def quick_alignment_distance(seq_i, seq_j, scoring_matrix=None, gap_bool=(True,True)):
	#
	# aln score
	#
	if scoring_matrix == None:
		aln = pairwise2.align.globalms(seq_i, seq_j, MATCH_NORMAL, XMATCH_NORMAL, GAP_OPEN, GAP_EXT, penalize_end_gaps=gap_bool)
		iden_score = MATCH_NORMAL*max(len(seq_i), len(seq_j))
	else:
		aln = pairwise2.align.globalds(seq_i, seq_j, scoring_matrix, GAP_OPEN, GAP_EXT, penalize_end_gaps=gap_bool)
		is1 = sum([scoring_matrix[(n,n)] for n in seq_i])
		is2 = sum([scoring_matrix[(n,n)] for n in seq_j])
		iden_score = max(is1, is2)
	aln_score = int(aln[0].score)
	#
	# rand score
	#
	rand_scores = []
	for k in range(RAND_SHUFFLE_COUNT):
		if scoring_matrix == None:
			rand_aln = pairwise2.align.globalms(shuffle_seq(seq_i), shuffle_seq(seq_j), MATCH_NORMAL, XMATCH_NORMAL, GAP_OPEN, GAP_EXT, penalize_end_gaps=gap_bool)
		else:
			rand_aln = pairwise2.align.globalds(shuffle_seq(seq_i), shuffle_seq(seq_j), scoring_matrix, GAP_OPEN, GAP_EXT, penalize_end_gaps=gap_bool)
		rand_scores.append(rand_aln[0].score)
	rand_score = int(np.mean(rand_scores))
	#
	# distance calculation
	#
	if rand_score >= aln_score:
		my_dist = MAX_SEQ_DIST
	else:
		my_dist = min(-np.log((aln_score - rand_score)/(iden_score - rand_score)), MAX_SEQ_DIST)
	my_dist = max(my_dist, MIN_SEQ_DIST)
	#
	# output
	#
	print(aln_score, rand_score, iden_score, '{0:0.3f}'.format(my_dist))
	return my_dist

#
#
#
def find_density_boundary(sequence, which_letter, win_size, dens_thresh, thresh_dir='below'):
	my_unknown = np.zeros(len(sequence))
	for j in range(len(sequence)):
		if sequence[j] == which_letter:
			my_unknown[j] = 1
	my_unknown_cum  = np.cumsum(my_unknown)
	my_unknown_dens = []
	first_pos_below_thresh = None
	pos_with_lowest_dens   = None	# use pos of min dense in the event that we never go below thresh
	if thresh_dir == 'below':
		lowest_dens_thus_far = 1.0
	elif thresh_dir == 'above':
		lowest_dens_thus_far = 0.0
	else:
		print('Error: find_density_boundary thresh_dir must be above or below')
		exit(1)
	for j in range(len(sequence) - win_size):
		my_unknown_dens.append(float(my_unknown_cum[j+win_size] - my_unknown_cum[j]) / win_size)
		if thresh_dir == 'below':
			if first_pos_below_thresh == None and my_unknown_dens[-1] <= dens_thresh:
				first_pos_below_thresh = j
			if my_unknown_dens[-1] < lowest_dens_thus_far:
				lowest_dens_thus_far = my_unknown_dens[-1]
				pos_with_lowest_dens = j
		else:
			if first_pos_below_thresh == None and my_unknown_dens[-1] >= dens_thresh:
				first_pos_below_thresh = j
			if my_unknown_dens[-1] > lowest_dens_thus_far:
				lowest_dens_thus_far = my_unknown_dens[-1]
				pos_with_lowest_dens = j
	####fig = mpl.figure(0)
	####mpl.plot(list(range(len(my_unknown_dens))), my_unknown_dens)
	####mpl.plot([first_pos_below_thresh, first_pos_below_thresh], [0,1], '--r')
	####mpl.plot([pos_with_lowest_dens, pos_with_lowest_dens], [0,1], '--b')
	####mpl.title(which_letter + ' ' + thresh_dir)
	####mpl.show()
	####mpl.close(fig)
	if first_pos_below_thresh == None:
		first_pos_below_thresh = pos_with_lowest_dens
	#
	return (sequence[:first_pos_below_thresh], sequence[first_pos_below_thresh:])

#
#
#
def find_cumulative_boundary(sequence, which_letters, cum_thresh=0.05, min_hits=100):
	hits = [1*(n in which_letters) for n in sequence]
	hits_cum = np.cumsum(hits)
	if hits_cum[-1] < min_hits:	# not enough hits to even bother trying
		return len(sequence)
	#
	hits_cum = hits_cum/hits_cum[-1]
	first_pos_below_thresh = len(sequence)
	for i in range(len(hits_cum)):
		if hits_cum[i] >= cum_thresh:
			first_pos_below_thresh = i
			break
	####fig = mpl.figure(0)
	####mpl.plot(list(range(len(hits_cum))), hits_cum)
	####mpl.plot([first_pos_below_thresh, first_pos_below_thresh], [0,1], '--r')
	####mpl.show()
	####mpl.close(fig)
	#
	return first_pos_below_thresh

#
#
#
def convert_colorvec_to_kmerhits(colorvecs, repeats_metadata):
	#
	[kmer_list, kmer_colors, kmer_letters, kmer_flags] = repeats_metadata
	#
	amino_2_kmer_ind = {}
	for i in range(len(kmer_colors)):
		amino_2_kmer_ind[kmer_letters[i]] = i
	#
	out_kmerhits = []
	for i in range(len(colorvecs)):
		current_block = colorvecs[i][0]
		current_start = 0
		out_kmerhits.append([[] for n in range(len(kmer_colors))])
		colorvecs[i] += UNKNOWN_LETTER
		for j in range(1,len(colorvecs[i])):
			if colorvecs[i][j] != current_block:
				if current_block != UNKNOWN_LETTER:
					my_ind = amino_2_kmer_ind[current_block]
					out_kmerhits[-1][my_ind].append((current_start, j))
				current_block = colorvecs[i][j]
				current_start = j
	return out_kmerhits

#
#
#
def denoise_colorvec(v, min_size=10, max_gap_fill=50, chars_to_delete=[], char_to_merge=''):
	blocks = []
	current_block = v[0]
	current_start = 0
	v += UNKNOWN_LETTER
	for i in range(1,len(v)):
		if v[i] != current_block:
			if current_block != UNKNOWN_LETTER:
				blocks.append((current_start, i, current_block))
			current_block = v[i]
			current_start = i
	#
	del_list = []
	for i in range(len(blocks)):
		if blocks[i][1] - blocks[i][0] < min_size and blocks[i][2] in chars_to_delete:
			del_list.append(i)
	del_list = sorted(del_list, reverse=True)
	for di in del_list:
		del blocks[di]
	#
	for i in range(len(blocks)-1,0,-1):
		if blocks[i][2] == blocks[i-1][2] and blocks[i][2] == char_to_merge:
			our_gap = blocks[i][0] - blocks[i-1][1]
			if our_gap <= max_gap_fill:
				blocks[i-1] = (blocks[i-1][0], blocks[i][1], blocks[i][2])
				del blocks[i]
	#
	v_out = [UNKNOWN_LETTER for n in v]
	for block in blocks:
		for i in range(block[0],block[1]):
			v_out[i] = block[2]
	return ''.join(v_out)

#
#	kmer_dat[i] = [[kmer1_hits, kmer2_hits, ...], tlen, tel-anchor-dist, read-orientation, readname, anchor_mapq]
#
#	repeats_metadata = [kmer_list, kmer_colors, kmer_letters, kmer_flags]
#
def cluster_tel_sequences(kmer_dat, repeats_metadata, my_chr, my_pos, dist_in=None, fig_name=None, msa_dir='', save_msa=None, train_prefix=None, tvr_truncate=3000, tree_cut=None, alignment_processes=8, muscle_exe='muscle', approx_rand=False):
	#
	[kmer_list, kmer_colors, kmer_letters, kmer_flags] = repeats_metadata
	#
	n_reads = len(kmer_dat)
	pq      = my_chr[-1]
	#
	if tree_cut == None:
		tree_cut = DEFAULT_TREECUT
	#
	canonical_letter = None
	denoise_chars    = []
	tvr_letters      = []
	for i in range(len(kmer_list)):
		if 'canonical' in kmer_flags[i]:
			canonical_letter = kmer_letters[i]
		if 'denoise' in kmer_flags[i]:
			denoise_chars.append(kmer_letters[i])
		if 'tvr' in kmer_flags[i]:
			tvr_letters.append(kmer_letters[i])
		if kmer_letters[i] == UNKNOWN_LETTER:
			print('Error: character A is reserved for unknown sequence')
			exit(1)
	if canonical_letter == None:
		print('Error: cluster_tel_sequences() received a kmer list that does not have CCCTAA or TTAGGG')
		exit(1)
	#
	# when generating consensus sequence for cluster: in ties, prioritize canonical, deprioritize unknown
	#
	char_score_adj = {canonical_letter:1, UNKNOWN_LETTER:-1}
	#
	# create color vector
	#
	all_colorvecs = []
	for i in range(n_reads):
		[my_kmer_hits, my_tlen, my_dbta, my_orr, my_rname, my_mapq] = kmer_dat[i]
		my_letters = [UNKNOWN_LETTER for n in range(my_tlen)]
		for ki in range(len(my_kmer_hits)):
			if len(my_kmer_hits[ki]):
				for kmer_span in my_kmer_hits[ki]:
					for j in range(kmer_span[0], kmer_span[1]):
						my_letters[j] = kmer_letters[ki]
		all_colorvecs.append(''.join(my_letters))
	#
	# identify subtel / tvr boundary
	#
	subtel_regions = []
	tvrtel_regions = []
	cleaned_colorvecs = []
	colorvecs_for_msa = []
	err_end_lens      = []
	for i in range(len(all_colorvecs)):
		#
		if pq == 'p':
			# identify subtel/tvr boundary based on density of unknown characters
			(seq_left, seq_right) = find_density_boundary(all_colorvecs[i][::-1], UNKNOWN_LETTER, UNKNOWN_WIN_SIZE, UNKNOWN_END_DENS, thresh_dir='below')
			# denoise tvr+tel section
			seq_right_denoise = denoise_colorvec(seq_right, chars_to_delete=denoise_chars, char_to_merge=canonical_letter)
			# remove ends of reads that might be sequencing artifacts, based on density of canonical characters
			(err_left, err_right) = find_density_boundary(seq_right_denoise[::-1], canonical_letter, CANON_WIN_SIZE, CANON_END_DENS, thresh_dir='above')
			#
			cleaned_colorvecs.append(err_right + seq_left[::-1])	# the entire denoised read (for plotting)
			err_end_lens.append(len(err_left))						# needed for adjusting offsets in plots
			#
			colorvecs_for_msa.append(seq_right_denoise[::-1])
			#
			subtel_regions.append(seq_left[::-1])					# subtel regions (currently not used for anything)
			tvrtel_regions.append(err_right)						# tvr sequence used for clustering
			if len(tvrtel_regions[-1]) > tvr_truncate:
				tvrtel_regions[-1] = tvrtel_regions[-1][-tvr_truncate:]
			#cleaned_colorvecs[-1] = tvrtel_regions[-1]
		#
		elif pq == 'q':
			# identify subtel/tvr boundary based on density of unknown characters
			(seq_left, seq_right) = find_density_boundary(all_colorvecs[i], UNKNOWN_LETTER, UNKNOWN_WIN_SIZE, UNKNOWN_END_DENS, thresh_dir='below')
			# denoise tvr+tel section
			seq_right_denoise = denoise_colorvec(seq_right, chars_to_delete=denoise_chars, char_to_merge=canonical_letter)
			# remove ends of reads that might be sequencing artifacts, based on density of canonical characters
			(err_left, err_right) = find_density_boundary(seq_right_denoise[::-1], canonical_letter, CANON_WIN_SIZE, CANON_END_DENS, thresh_dir='above')
			#
			cleaned_colorvecs.append(seq_left + err_right[::-1])	# the entire denoised read (for plotting)
			err_end_lens.append(len(err_left))						# needed for adjusting offsets in plots
			#
			colorvecs_for_msa.append(seq_right_denoise)
			#
			subtel_regions.append(seq_left)							# subtel regions (currently not used for anything)
			tvrtel_regions.append(err_right[::-1])					# tvr sequence used for clustering
			if len(tvrtel_regions[-1]) > tvr_truncate:
				tvrtel_regions[-1] = tvrtel_regions[-1][:tvr_truncate]
			#cleaned_colorvecs[-1] = tvrtel_regions[-1]
	#
	# trivial case
	#
	if n_reads == 1:
		return [ [[0]], [[kmer_dat[0][5]]], [[0]], tvrtel_regions, cleaned_colorvecs, err_end_lens ]
	#
	# PAIRWISE ALIGNMENT OF ALL SEQUENCES
	#
	if dist_in == None or exists_and_is_nonzero(dist_in) == False:
		all_indices = [[] for n in range(alignment_processes)]
		k = 0
		for i in range(n_reads):
			for j in range(i+1,n_reads):
				all_indices[k%alignment_processes].append((i,j))
				k += 1
		#
		# align tvr + tel regions
		#
		manager     = multiprocessing.Manager()
		tvrtel_dist = manager.dict()
		processes   = []
		train_handl = []
		for i in range(alignment_processes):
			if train_prefix == None or exists_and_is_nonzero(train_prefix+'_tvrtel.tsv'):
				train_handl.append(None)
			else:
				job_fn = train_prefix + '.' + str(i+1).zfill(len(str(alignment_processes))) + '_tvrtel.tsv'
				train_handl.append(job_fn)
				f = open(job_fn, 'w')
				f.close()
			p = multiprocessing.Process(target=parallel_alignment_job, args=(all_indices[i], tvrtel_regions, (True,True), pq, tvrtel_dist, None, train_handl[-1], approx_rand))
			processes.append(p)
		for i in range(alignment_processes):
			processes[i].start()
		for i in range(alignment_processes):
			processes[i].join()
		if train_prefix != None and exists_and_is_nonzero(train_prefix+'_tvrtel.tsv') == False:
			os.system('cat ' + ' '.join(train_handl) + ' > ' + train_prefix + '_tvrtel.tsv')
			for n in train_handl:
				os.system('rm ' + n)
		#
		# combine distances for final distance matrix
		#
		dist_matrix = np.zeros((n_reads,n_reads))
		for (i,j) in tvrtel_dist.keys():
			ij_dist = tvrtel_dist[(i,j)]
			dist_matrix[i,j] = ij_dist
			dist_matrix[j,i] = ij_dist
		dist_norm    = max(np.max(dist_matrix), MIN_MSD)
		dist_matrix /= dist_norm
		#
		if dist_in != None:
			np.save(dist_in, dist_matrix)
		#
	else:
		dist_matrix = np.load(dist_in, allow_pickle=True)
	#
	# hierarchal clustering + dendrogram plotting
	#
	dist_array = squareform(dist_matrix)
	Zread      = linkage(dist_array, method='ward')
	#
	if fig_name != None:
		fig = mpl.figure(3, figsize=(8,6))
		mpl.rcParams.update({'font.size': 16, 'font.weight':'bold', 'lines.linewidth':2.0})
		dendrogram(Zread, color_threshold=tree_cut)
		mpl.axhline(y=[tree_cut], linestyle='dashed', color='black', alpha=0.7)
		mpl.xlabel('read #')
		mpl.ylabel('distance')
		mpl.title(my_chr + ' : ' + str(my_pos))
		mpl.tight_layout()
		mpl.savefig(fig_name)
		mpl.close(fig)
	#
	assignments = fcluster(Zread, tree_cut, 'distance').tolist()
	by_class = {}
	for i in range(len(assignments)):
		if assignments[i] not in by_class:
			by_class[assignments[i]] = []
		by_class[assignments[i]].append(i)
	out_clust = sorted([(len(by_class[k]), sorted(by_class[k])) for k in by_class.keys()], reverse=True)
	out_clust = [n[1] for n in out_clust]
	for i in range(len(out_clust)):	# sort by length
		out_clust[i] = sorted([(kmer_dat[n][1],n) for n in out_clust[i]], reverse=True)
		out_clust[i] = [n[1] for n in out_clust[i]]
	#
	# do MSA to get a consensus sequence
	#
	out_consensus = []
	if save_msa == None or exists_and_is_nonzero(save_msa) == False:
		for i in range(len(out_clust)):
			max_seq_len = max([len(colorvecs_for_msa[n]) for n in out_clust[i]])
			for ci in out_clust[i]:
				buff_seq = canonical_letter*(max_seq_len - len(colorvecs_for_msa[ci]))
				if pq == 'p':
					colorvecs_for_msa[ci] = buff_seq + colorvecs_for_msa[ci]
				elif pq == 'q':
					colorvecs_for_msa[ci] = colorvecs_for_msa[ci] + buff_seq
		# 
		for i in range(len(out_clust)):
			if len(out_clust[i]) == 1:
				out_consensus.append(colorvecs_for_msa[out_clust[i][0]])
			else:
				clust_seq = [colorvecs_for_msa[n] for n in out_clust[i]]
				[msa_seq, consensus_seq] = get_muscle_msa(clust_seq, muscle_exe, working_dir=msa_dir, char_score_adj=char_score_adj)
				out_consensus.append(consensus_seq)
		#
		if save_msa != None:
			f = open(save_msa,'w')
			for i in range(len(out_consensus)):
				f.write('>msa-' + str(i+1).zfill(2) + '\n')
				f.write(out_consensus[i] + '\n')
			f.close()
	else:
		my_reader = TG_Reader(save_msa, verbose=False)
		while True:
			read_dat = my_reader.get_next_read()
			if not read_dat[0]:
				break
			out_consensus.append(read_dat[1])
		my_reader.close()

	#
	# prune subtel from consensus
	#
	for i in range(len(out_consensus)):
		if pq == 'p':
			(cons_left, cons_right) = find_density_boundary(out_consensus[i][::-1], UNKNOWN_LETTER, UNKNOWN_WIN_SIZE, UNKNOWN_END_DENS, thresh_dir='below')
			out_consensus[i] = cons_right[::-1]
		elif pq == 'q':
			(cons_left, cons_right) = find_density_boundary(out_consensus[i], UNKNOWN_LETTER, UNKNOWN_WIN_SIZE, UNKNOWN_END_DENS, thresh_dir='below')
			out_consensus[i] = cons_right
	#
	# identify tvr/tel boundary from consensus sequences
	#
	out_tvr_tel_boundaries = []
	for i in range(len(out_consensus)):
		if len(out_consensus[i]) > MAX_TVR_LEN:
			if pq == 'p':
				current_cons = canonical_letter*(len(out_consensus[i])-MAX_TVR_LEN) + out_consensus[i][-MAX_TVR_LEN:]
			elif pq == 'q':
				current_cons = out_consensus[i][:MAX_TVR_LEN] + canonical_letter*(len(out_consensus[i])-MAX_TVR_LEN)
		else:
			current_cons = out_consensus[i]
		denoised_consensus = denoise_colorvec(current_cons, chars_to_delete=tvr_letters, min_size=TVR_CANON_FILT_PARAMS_STRICT[0], char_to_merge=canonical_letter)
		if pq == 'q':
			denoised_consensus = denoised_consensus[::-1]
		tel_boundary = find_cumulative_boundary(denoised_consensus, tvr_letters, cum_thresh=TVR_CANON_FILT_PARAMS_STRICT[1], min_hits=TVR_CANON_FILT_PARAMS_STRICT[2])
		#
		#
		#
		if tel_boundary == len(out_consensus[i])+1:	# failed to find tel boundary, try again with lenient params
			if len(out_consensus[i]) > MAX_TVR_LEN_SHORT:
				if pq == 'p':
					current_cons = canonical_letter*(len(out_consensus[i])-MAX_TVR_LEN_SHORT) + out_consensus[i][-MAX_TVR_LEN_SHORT:]
				elif pq == 'q':
					current_cons = out_consensus[i][:MAX_TVR_LEN_SHORT] + canonical_letter*(len(out_consensus[i])-MAX_TVR_LEN_SHORT)
			else:
				current_cons = out_consensus[i]
			denoised_consensus = denoise_colorvec(current_cons, chars_to_delete=tvr_letters, min_size=TVR_CANON_FILT_PARAMS_LENIENT[0], char_to_merge=canonical_letter)
			if pq == 'q':
				denoised_consensus = denoised_consensus[::-1]
			tel_boundary = find_cumulative_boundary(denoised_consensus, tvr_letters, cum_thresh=TVR_CANON_FILT_PARAMS_LENIENT[1], min_hits=TVR_CANON_FILT_PARAMS_LENIENT[2])
			#print('LENIENT TVR BOUNDARY (cluster '+str(i)+'):', len(out_consensus[i]) - tel_boundary + 1)
		out_tvr_tel_boundaries.append(len(out_consensus[i]) - tel_boundary + 1)
	#
	# lets use tvr/subtel boundary as offset instead of msa offset
	#
	all_subtel_lens = [len(subtel_regions[n]) for n in range(len(subtel_regions))]
	longest_subtel  = max(all_subtel_lens)
	out_adj         = []
	for i in range(len(out_clust)):
		my_subtel_lens    = [len(subtel_regions[n]) for n in out_clust[i]]
		longest_subtel_cl = max(my_subtel_lens)
		clust_subtel_adj  = longest_subtel - longest_subtel_cl
		out_adj.append([longest_subtel_cl - n + clust_subtel_adj for n in my_subtel_lens])
	#
	# output mapping quality as well (not sure what I'll do with this yet, maybe filtering?)
	#
	out_mapq = []
	for n in out_clust:
		out_mapq.append([kmer_dat[m][5] for m in n])
	#
	print(out_clust)
	print(all_subtel_lens)
	print(out_tvr_tel_boundaries)
	#
	return [ out_clust, out_mapq, out_adj, all_subtel_lens, out_consensus, cleaned_colorvecs, err_end_lens, out_tvr_tel_boundaries ]

#
#
#
def cluster_consensus_tvr(sequences, repeats_metadata, dist_in=None, fig_name=None, samp_labels=None, aln_mode='ms', linkage_method='complete', tree_cut=4.30, alignment_processes=8, job=(1,1), dendrogram_height=12):
	#
	n_seq = len(sequences)
	#
	[kmer_list, kmer_colors, kmer_letters, kmer_flags] = repeats_metadata
	#
	canonical_letter = None
	for i in range(len(kmer_list)):
		if 'canonical' in kmer_flags[i]:
			canonical_letter = kmer_letters[i]
		if kmer_letters[i] == UNKNOWN_LETTER:
			print('Error: character A is reserved for unknown sequence')
			exit(1)
	if canonical_letter == None:
		print('Error: cluster_consensus_tvr() received a kmer list that does not have CCCTAA or TTAGGG')
		exit(1)
	#
	if dist_in == None or exists_and_is_nonzero(dist_in) == False:
		dist_matrix = np.zeros((n_seq,n_seq))
		pw2_gap     = (False, False)
		all_indices = [[] for n in range(alignment_processes)]
		k = 0
		for i in range(n_seq):
			for j in range(i+1,n_seq):
				all_indices[k%alignment_processes].append((i,j))
				k += 1
		#
		#	scoring matrix
		#
		letters = AMINO
		scoring_matrix = {}
		for i in range(len(letters)):
			for j in range(len(letters)):
				if i == j:
					scoring_matrix[(letters[i],letters[j])] = MATCH_NORMAL
				else:
					scoring_matrix[(letters[i],letters[j])] = XMATCH_NORMAL
					scoring_matrix[(letters[j],letters[i])] = XMATCH_NORMAL
		for i in range(len(letters)):
			scoring_matrix[(letters[i],canonical_letter)] = XMATCH_CANON
			scoring_matrix[(canonical_letter,letters[i])] = XMATCH_CANON
		for i in range(len(letters)):
			scoring_matrix[(letters[i],UNKNOWN_LETTER)] = XMATCH_UNKNOWN
			scoring_matrix[(UNKNOWN_LETTER,letters[i])] = XMATCH_UNKNOWN
		scoring_matrix[(canonical_letter, canonical_letter)] = MATCH_CANON		# reduced award for matching canonical
		scoring_matrix[(UNKNOWN_LETTER, UNKNOWN_LETTER)]     = MATCH_UNKNOWN	# no reward for matching unknown regions
		#
		#	even more parallelization! Any problem can be solved by throwing tons of CPU at it.
		#
		if job[1] > 1:
			my_job = job[0]-1
			chunks = job[1]
			for i in range(alignment_processes):
				chunksize = int(len(all_indices[i])/chunks)
				chunks_by_job = []
				for j in range(chunks):
					if j == chunks-1:
						chunks_by_job.append(all_indices[i][j*chunksize:])
					else:
						chunks_by_job.append(all_indices[i][j*chunksize:(j+1)*chunksize])
				all_indices[i] = [n for n in chunks_by_job[my_job]]
		#
		manager     = multiprocessing.Manager()
		return_dict = manager.dict()
		processes   = []
		for i in range(alignment_processes):
			if aln_mode == 'ms':
				p = multiprocessing.Process(target=parallel_alignment_job, args=(all_indices[i], sequences, pw2_gap, 'q', return_dict, None, None, True))
			elif aln_mode == 'ds':
				p = multiprocessing.Process(target=parallel_alignment_job, args=(all_indices[i], sequences, pw2_gap, 'q', return_dict, scoring_matrix, None, False))
			processes.append(p)
		for i in range(alignment_processes):
			processes[i].start()
		for i in range(alignment_processes):
			processes[i].join()
		#
		dist_matrix = np.zeros((n_seq,n_seq))
		for (i,j) in return_dict.keys():
			dist_matrix[i,j] = return_dict[(i,j)]
			dist_matrix[j,i] = return_dict[(i,j)]
		if job[1] > 1:
			partial_dist_fn = dist_in[:-4] + '_job' + str(job[0]).zfill(3) + '.npy'
			np.save(partial_dist_fn, dist_matrix)
		else:
			dist_norm    = max(np.max(dist_matrix), MIN_MSD)
			dist_matrix /= dist_norm
			if dist_in != None:
				np.save(dist_in, dist_matrix)
	else:
		dist_matrix = np.load(dist_in, allow_pickle=True)
	#
	if job[1] == 1 or job[0] == 0:
		d_arr = squareform(dist_matrix)
		Zread = linkage(d_arr, method=linkage_method)
		#
		if fig_name != None:
			mpl.rcParams.update({'font.size': 16, 'font.weight':'bold'})
			#
			fig = mpl.figure(3, figsize=(8,dendrogram_height))
			dendro_dat = dendrogram(Zread, orientation='left', labels=samp_labels, color_threshold=tree_cut)
			mpl.axvline(x=[tree_cut], linestyle='dashed', color='black', alpha=0.7)
			mpl.xlabel('distance')
			#
			mpl.tight_layout()
			mpl.savefig(fig_name, dpi=200)	# default figure dpi = 100
			mpl.close(fig)
		#
		labels_fromtop = dendro_dat['ivl'][::-1]
		#
		assignments = fcluster(Zread, tree_cut, 'distance').tolist()
		by_class = {}
		for i in range(len(assignments)):
			if assignments[i] not in by_class:
				by_class[assignments[i]] = []
			by_class[assignments[i]].append(i)
		out_clust = sorted([(len(by_class[k]), sorted(by_class[k])) for k in by_class.keys()], reverse=True)
		out_clust = [n[1] for n in out_clust]
		#
		return out_clust
	#
	return None
