import re
import os
import copy

import numpy as np

LEXICO_2_IND = {'chr1':1, 'chr2':2, 'chr3':3, 'chr10':10, 'chr11':11, 'chr12':12, 'chr19':19, 'chr20':20,
				'chr4':4, 'chr5':5, 'chr6':6, 'chr13':13, 'chr14':14, 'chr15':15, 'chr21':21, 'chr22':22,
				'chr7':7, 'chr8':8, 'chr9':9, 'chr16':16, 'chr17':17, 'chr18':18, 'chrX' :23, 'chrY' :24, 'chrM' :25}

T2T_CHROMSIZE = {'chr1':248387328,
                 'chr2':242696752,
                 'chr3':201105948,
                 'chr4':193574945,
                 'chr5':182045439,
                 'chr6':172126628,
                 'chr7':160567428,
                 'chr8':146259331,
                 'chr9':150617247,
                 'chr10':134758134,
                 'chr11':135127769,
                 'chr12':133324548,
                 'chr13':113566686,
                 'chr14':101161492,
                 'chr15':99753195,
                 'chr16':96330374,
                 'chr17':84276897,
                 'chr18':80542538,
                 'chr19':61707364,
                 'chr20':66210255,
                 'chr21':45090682,
                 'chr22':51324926,
                 'chrX':154259566,
                 'chrY':62460029}

LARGE_NUMBER = int(1e9)

def exists_and_is_nonzero(fn):
	if os.path.isfile(fn):
		if os.path.getsize(fn) > 0:
			return True
	return False

def makedir(d):
	if not os.path.isdir(d):
		os.mkdir(d)

def rm(fn):
	if os.path.isdir(fn):
		os.rmdir(fn)
	elif os.path.isfile(fn):
		os.remove(fn)

#
#
RC_DICT = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
def RC(s):
	return ''.join([RC_DICT[n] for n in s[::-1]])

#
#
REF_CHAR  = 'MX=D'
READ_CHAR = 'MX=I'
CLIP_CHAR = 'SH'
def parse_cigar(cigar):
	letters = re.split(r"\d+",cigar)[1:]
	numbers = [int(n) for n in re.findall(r"\d+",cigar)]
	startPos = 0
	if letters[0] in CLIP_CHAR:
		startPos = numbers[0]
	endClip = 0
	if len(letters) > 1 and letters[-1] in CLIP_CHAR:
		endClip = numbers[-1]
	adj  = 0
	radj = 0
	for i in range(len(letters)):
		if letters[i] in REF_CHAR:
			adj += numbers[i]
		if letters[i] in READ_CHAR:
			radj += numbers[i]
	return (startPos, adj, radj, endClip)

#
#
#
def parse_read(splt):
	#
	cigar = splt[5]
	flag  = int(splt[1])
	ref   = splt[2]
	pos   = int(splt[3])
	rdat  = splt[9]
	rnm   = splt[0]
	mapq  = int(splt[4])
	#
	is_unmapped = False
	if ref == '*' or flag&4:
		is_unmapped = True
	#
	orientation = 'FWD'
	if flag&16:
		orientation = 'REV'
	#
	if ref[:3] == 'chr':	# assume format: "chr1_p"
		ref_key = ref.split('_')[0]
		if ref_key in LEXICO_2_IND:
			ref_key = LEXICO_2_IND[ref_key]
		else:
			print('unknown ref name:', ref, '-->', ref_key)
			exit(1)
	#
	elif ref[:3] == 'alt':	# assume format: "alt1_p_0" OR "alt_1p_0" because I'm sloppy
		ref_key = ''.join(ref.split('_')[:-1])
		ref_key = ref_key.replace('alt', 'chr')[:-1]
		alt_mul = int(ref.split('_')[-1]) + 1
		if ref_key in LEXICO_2_IND:
			ref_key = LEXICO_2_IND[ref_key] + alt_mul*len(LEXICO_2_IND)
		else:
			print('unknown ref name:', ref, '-->', ref_key)
			exit(1)
	#
	elif ref[:3] == 'tel':
		ref_key = LARGE_NUMBER
	#
	elif is_unmapped == True:
		ref_key = LARGE_NUMBER + 1
	#
	else:
		print('unknown ref name:', ref)
		exit(1)
	#
	(read_pos_1, read_pos_2) = (None, None)
	(pos1, pos2) = (None, None)
	if is_unmapped == False:
		cig_dat    = parse_cigar(cigar)
		read_pos_1 = cig_dat[0]
		read_pos_2 = cig_dat[0] + cig_dat[2]
		read_len   = cig_dat[0] + cig_dat[2] + cig_dat[3]
		#
		if orientation == 'REV':
			[read_pos_1, read_pos_2] = [read_len - read_pos_2, read_len - read_pos_1]
		if orientation == 'FWD':
			pos1, pos2 = pos, pos + cig_dat[1]
		elif orientation == 'REV':
			pos1, pos2 = pos + cig_dat[1], pos
	#
	return [rnm, ref_key, pos, read_pos_1, read_pos_2, ref, pos1, pos2, orientation, mapq, rdat]

#
# trim repeated matches in the same manner as pbmm2 (only affects read_pos coords)
#
def repeated_matches_trimming(alns, min_read_span_after_trimming=200, strategy='mapq', print_debug=False):
	if print_debug:
		print('- matches trimming')
		for n in alns:
			print(n[:7], 'rdat len:', len(n[7]))
		print()
	r_coords = [[alns[n][0], alns[n][1], n] for n in range(len(alns)) if (alns[n][0] != None and alns[n][1] != None)]
	# don't try to trim unmapped reads
	if len(r_coords) == 0:
		return alns
	clust    = cluster_ranges(r_coords)
	any_lap  = any([len(n) > 1 for n in clust])
	# we have nothing to trim
	if any_lap == False:
		return alns
	#
	while any_lap:
		flat_clust = []
		#
		for i in range(len(clust)):
			if len(clust[i]) > 1:
				max_span = (clust[i][0][0], clust[i][0][1], alns[clust[i][0][2]][6])
				for j in range(1,len(clust[i])):
					if strategy == 'largest' and clust[i][j][1] - clust[i][j][0] > max_span[1] - max_span[0]:
						max_span = (clust[i][j][0], clust[i][j][1], alns[clust[i][j][2]][6])
					elif strategy == 'mapq' and alns[clust[i][j][2]][6] > max_span[2]:
						max_span = (clust[i][j][0], clust[i][j][1], alns[clust[i][j][2]][6])
				if print_debug:
					print('MAX_SPAN', max_span)
				#
				max_found = False
				del_list  = []
				for j in range(len(clust[i])):
					(x, y, q) = (clust[i][j][0], clust[i][j][1], alns[clust[i][j][2]][6])
					# we are max (use first, if multiple)
					if x == max_span[0] and y == max_span[1]:
						if strategy == 'largest':
							if max_found:
								del_list.append(j)
							else:
								max_found = True
						elif strategy == 'mapq':
							if q == max_span[2] and max_found == False:
								max_found = True
							else:
								del_list.append(j)
					# are we completely consumed by max?
					elif x >= max_span[0] and y <= max_span[1]:
						del_list.append(j)
					# do we overhang on the left?
					elif x < max_span[0] and y >= max_span[0]:
						clust[i][j][1] = max_span[0]-1
					# do we overhang on the right?
					elif x <= max_span[1] and y > max_span[1]:
						clust[i][j][0] = max_span[1]+1
					# did the trimming make us too short?
					if clust[i][j][1] - clust[i][j][0] < min_read_span_after_trimming:
						del_list.append(j)
				#
				del_list = sorted(list(set(del_list)), reverse=True)
				for di in del_list:
					del clust[i][di]
			flat_clust.extend(clust[i])
		#
		if len(flat_clust) == 0:
			return []
		#
		flat_clust = sorted(flat_clust)
		clust      = cluster_ranges(flat_clust)
		any_lap    = any([len(n) > 1 for n in clust])
	#
	alns_out = copy.deepcopy(alns)
	del_list = []
	keep_ind = {n[2]:(n[0],n[1]) for n in flat_clust}
	for i in range(len(alns_out)):
		if i in keep_ind:
			alns_out[i][0] = keep_ind[i][0]
			alns_out[i][1] = keep_ind[i][1]
		else:
			del_list.append(i)
	del_list = sorted(list(set(del_list)), reverse=True)
	for di in del_list:
		del alns_out[di]
	#
	if print_debug:
		for n in alns_out:
			print(n[:7])
		print()
	return alns_out

#
# cluster a sorted list
#
# "which_val = None" - assume input is list of values to directly sort
# "which_val = 1"    - assume input is list of tuples, use index 1 to sort
#
def cluster_list(l, delta, which_val=None):
	if which_val == None:
		prev_val = l[0]
	else:
		prev_val = l[0][which_val]
	out_list    = [[l[0]]]
	current_ind = 0
	for n in l[1:]:
		if which_val == None:
			my_dist = n - prev_val
		else:
			my_dist = n[which_val] - prev_val
		if my_dist <= delta:
			out_list[current_ind].append(n)
		else:
			current_ind += 1
			out_list.append([])
			out_list[current_ind].append(n)
		if which_val == None:
			prev_val = n
		else:
			prev_val = n[which_val]
	return out_list

#
#
#
def cluster_ranges(l):
	c = [[l[0]]]
	for i in range(1,len(l)):
		found_a_home = False
		for j in range(len(c)):
			for k in range(len(c[j])):
				if l[i][0] <= c[j][k][1] and l[i][1] >= c[j][k][0]:
					c[j].append(l[i])
					found_a_home = True
					break
			if found_a_home:
				break
		if found_a_home == False:
			c.append([l[i]])
	return c

#
#
#
def posmax(seq):
	m = seq[0]
	index = 0
	for i,x in enumerate(seq):
		if x > m:
			m = x
			index = i
	return index
