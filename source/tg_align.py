import numpy as np
import time

from Bio.Align          import PairwiseAligner, substitution_matrices
from collections        import Counter, defaultdict
from concurrent.futures import as_completed, FIRST_COMPLETED, ProcessPoolExecutor, wait
from itertools          import combinations, zip_longest

from source.tg_util import shuffle_seq


# we're going to pretend each kmer color is an amino acid, so alignment tools behave themselves
AMINO = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
AMINO_2_IND = {AMINO[i]:i for i in range(len(AMINO))}
UNKNOWN = 'A'
GAP_CHR = '-'

# when comparing sequences of different lengths, choose min(len(seq1), len(seq2)), but only go this low
MIN_VIABLE_SEQ_LEN = 1000
# the log-based distance function can go to infinity so lets set an upper bound on it
MAX_SEQ_DIST = 10.0
# similarly, lets choose a small number as the minimum to prevent numerical weirdness from giving us negative values
MIN_SEQ_DIST = 0.0001
# if aln_score is < - IDEN_SCORE_SCALAR_CHEAT * iden_score, then assign MAX_SEQ_DIST
IDEN_SCORE_SCALAR_CHEAT = 0.900

# use this many reads when creating MSAs for consensus sequences (chooses longest reads)
MAX_MSA_READCOUNT = 10


def get_scoring_matrix(canon_char, which_type='tvr'):
    m = substitution_matrices.Array(''.join(AMINO) + GAP_CHR, dims=2)
    for c1 in AMINO:
        m[c1,c1]                                      = 5.  # normal match
        for c2 in AMINO:
            if c1 != c2:
                m[c1,c2]                              = -4. # normal mismatch
    m[UNKNOWN,UNKNOWN]                                = 2.  # matching unknown
    m[GAP_CHR,GAP_CHR]                                = 0.  # matching gap
    for c1 in AMINO:
        m[c1,canon_char] = m[canon_char,c1]           = -4. # canon mismatch
        m[c1,UNKNOWN] = m[UNKNOWN,c1]                 = -2. # unknown mismatch
        m[c1,GAP_CHR] = m[GAP_CHR,c1]                 = -4. # gap mismatch
    #
    if which_type == 'tvr':
        m[canon_char,canon_char]                      = 0.  # matching canonical
        m[canon_char,UNKNOWN] = m[UNKNOWN,canon_char] = -2. # canon = unknown
        m[UNKNOWN,GAP_CHR] = m[GAP_CHR,UNKNOWN]       = -2. # gap = unknown
    elif which_type == 'consensus':
        m[canon_char,canon_char]                      = 1.  # matching canonical
        m[canon_char,UNKNOWN] = m[UNKNOWN,canon_char] = -2. # canon = unknown
        m[UNKNOWN,GAP_CHR] = m[GAP_CHR,UNKNOWN]       = -2. # gap = unknown
    elif which_type == 'msa':
        m[canon_char,canon_char]                      = 0.  # matching canonical
        m[canon_char,UNKNOWN] = m[UNKNOWN,canon_char] = 1.  # canon = unknown
        m[UNKNOWN,GAP_CHR] = m[GAP_CHR,UNKNOWN]       = -2. # gap = unknown
    elif which_type == 'msa_refinement':
        for c1 in AMINO:
            m[c1,GAP_CHR] = m[GAP_CHR,c1]             = 0.  # gap mismatch
        m[canon_char,canon_char]                      = 2.  # matching canonical
        m[canon_char,UNKNOWN] = m[UNKNOWN,canon_char] = 1.  # canon = unknown
        m[UNKNOWN,GAP_CHR] = m[GAP_CHR,UNKNOWN]       = 1.  # gap = unknown
    else:
        print('Error: unknown which_type')
        exit(1)
    #
    return m


def get_aligner_object(scoring_matrix=None, gap_bool=(True,True), which_type='tvr'):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.alphabet = ''.join(AMINO) + GAP_CHR
    if which_type == 'tvr':
        gap_open = -4.
        gap_extend = -4.
    elif which_type == 'msa':
        gap_open = -16.
        gap_extend = -4.
    else:
        print('Error: unknown which_type')
        exit(1)
    aligner.match_score = 5.
    aligner.mismatch_score = -4.
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    #
    if scoring_matrix is not None:
        aligner.substitution_matrix = scoring_matrix
    if gap_bool[0]:
        aligner.left_open_gap_score   = gap_open
        aligner.left_extend_gap_score = gap_extend
    else:
        aligner.left_open_gap_score   = 0.0
        aligner.left_extend_gap_score = 0.0
    if gap_bool[1]:
        aligner.right_open_gap_score   = gap_open
        aligner.right_extend_gap_score = gap_extend
    else:
        aligner.right_open_gap_score   = 0.0
        aligner.right_extend_gap_score = 0.0
    return aligner


def tvr_distance(tvr_i, tvr_j, aligner, adjust_lens=True, min_viable=True, randshuffle=3):
    if adjust_lens:
        min_len = min(len(tvr_i), len(tvr_j))
        if min_viable:
            min_len = max(min_len, MIN_VIABLE_SEQ_LEN)
        seq_i = tvr_i[:min_len]
        seq_j = tvr_j[:min_len]
    else:
        seq_i = tvr_i
        seq_j = tvr_j
    # aln score
    aln_score = aligner.score(seq_i, seq_j)
    # iden score
    c1 = Counter(seq_i)
    c2 = Counter(seq_j)
    if aligner.substitution_matrix is None:
        iden_score = aligner.match_score * ((len(seq_i)+len(seq_j)) / 2.)
    else:
        is1 = sum([c1[n]*aligner.substitution_matrix[n,n] for n in c1.keys()])
        is2 = sum([c2[n]*aligner.substitution_matrix[n,n] for n in c2.keys()])
        iden_score = (is1+is2) / 2.
    # a cheat for performance reasons:
    # - if our aln_score is very low that almost always results in a MAX_SEQ_DIST score
    if aln_score < -(IDEN_SCORE_SCALAR_CHEAT * iden_score):
        return MAX_SEQ_DIST
    # rand score
    rand_scores = []
    for k in range(randshuffle):
        rand_scores.append(aligner.score(shuffle_seq(seq_i), shuffle_seq(seq_j)))
    rand_score_shuffle = np.mean(rand_scores)
    # convert to distance, bounds checking
    if rand_score_shuffle >= aln_score:
        dist_out = MAX_SEQ_DIST
    else:
        dist_out = min(-np.log((aln_score - rand_score_shuffle) / (iden_score - rand_score_shuffle)), MAX_SEQ_DIST)
    dist_out = max(dist_out, MIN_SEQ_DIST)
    #print(aln_score, iden_score, rand_score_shuffle, dist_out)
    return dist_out


def get_dist_matrix_parallel(sequences, aligner, adjust_lens, min_viable, rand_shuffle_count, max_workers=4, max_pending=100, tasks_per_worker=5000, print_progress=False):
    n_seq = len(sequences)
    tasks = combinations(range(n_seq),2)
    dist_matrix = np.zeros((n_seq,n_seq), dtype='<f4')
    tt = time.perf_counter()
    current_i = -1
    while True:
        # create a new ProcessPoolExecutor periodically to force process recycling
        # -- this is inelegant, but needed because of a memory leak in PairwiseAligner
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            pending_futures = {}
            tasks_completed = 0
            while tasks_completed < tasks_per_worker * max_workers:
                while len(pending_futures) < max_pending:
                    try:
                        (i,j) = next(tasks)
                        future = executor.submit(tvr_distance, sequences[i], sequences[j], aligner, adjust_lens, min_viable, rand_shuffle_count)
                        pending_futures[future] = (i,j)
                    except StopIteration:
                        break
                if not pending_futures:
                    return dist_matrix
                #
                done, _ = wait(pending_futures, return_when=FIRST_COMPLETED)
                for future in done:
                    my_tvr_dist = future.result()
                    (i,j) = pending_futures.pop(future)
                    dist_matrix[i,j] = dist_matrix[j,i] = my_tvr_dist
                    tasks_completed += 1
                    #
                    if i > current_i:
                        current_i = i
                        if print_progress:
                            print(f'progress: {i+1}/{n_seq} ({int(time.perf_counter() - tt)} sec)')


def remove_gap_columns(alignment):
    aln_out = [col for col in zip(*alignment) if set(col) != {'-'}]
    aln_out = [''.join(row) for row in zip(*aln_out)]
    return aln_out


def pairwise_align(seq1, seq2, aligner):
    alignment = aligner.align(seq1, seq2)[0]
    return str(alignment[0]), str(alignment[1])


def progressive_alignment(sequences, distance_matrix, aligner, num_processes=1):

    def align_to_profile(sequence, profile):
        alignments = [pairwise_align(seq, sequence, aligner) for seq, _ in profile]
        max_len = max(len(n[0]) for n in alignments)
        padded_profile = [(a[0].ljust(max_len, '-'), idx) for (a, idx) in zip(alignments, [idx for _, idx in profile])]
        padded_profile.append((alignments[0][1].ljust(max_len, '-'), len(sequences) - 1))
        return padded_profile

    def align_to_profile_parallel(sequence, profile):
        with ProcessPoolExecutor(max_workers=num_processes) as executor:
            future_alignments = [executor.submit(pairwise_align, seq, sequence, aligner) for seq, _ in profile]
            alignments = [future.result() for future in as_completed(future_alignments)]
        max_len = max(len(n[0]) for n in alignments)
        padded_profile = [(a[0].ljust(max_len, '-'), idx) for (a, idx) in zip(alignments, [idx for _, idx in profile])]
        padded_profile.append((alignments[0][1].ljust(max_len, '-'), len(sequences) - 1))
        return padded_profile

    n_sequences = len(sequences)
    if n_sequences <= 0:
        print('Error: progressive_alignment() received empty alignment')
        exit(1)
    # we need at least 3 sequences to do anything meaningful
    if n_sequences < 3:
        return sequences
    # find the two closest sequences
    min_distance = float('inf')
    closest_pair = None
    for i in range(n_sequences):
        for j in range(i+1, n_sequences):
            if distance_matrix[i][j] < min_distance:
                min_distance = distance_matrix[i][j]
                closest_pair = (i, j)
    seq1, seq2 = pairwise_align(sequences[closest_pair[0]], sequences[closest_pair[1]], aligner)
    alignment = [(seq1, closest_pair[0]), (seq2, closest_pair[1])]
    remaining = [(seq, i) for i, seq in enumerate(sequences) if i not in closest_pair]
    #
    while remaining:
        # get closest sequence to the current alignment
        distances = [np.mean([distance_matrix[i][j] for _, j in alignment]) for _, i in remaining]
        closest_index = np.argmin(distances)
        closest_seq, original_index = remaining.pop(closest_index)
        # align the closest sequence to the current profile
        if num_processes > 1:
            alignment = align_to_profile_parallel(closest_seq, alignment)
        else:
            alignment = align_to_profile(closest_seq, alignment)
    # sort based on the original sequence indices
    alignment.sort(key=lambda x: x[1])
    return remove_gap_columns([seq for seq, _ in alignment])


def iterative_refinement(initial_alignment, aligner, max_iterations=10, improvement_threshold=0.01):

    def create_profile(alignment):
        profile = {}
        for i in range(len(alignment[0])):
            column = [seq[i] for seq in alignment]
            profile[i] = {char: column.count(char) / len(column) for char in set(column)}
        return profile

    def score_against_profile(sequence, profile):
        score = 0
        for i, char in enumerate(sequence):
            if i in profile:
                for prof_char, freq in profile[i].items():
                    if char != '-' and prof_char != '-':
                        if aligner.substitution_matrix is None:
                            if char == prof_char:
                                score += aligner.match_score
                            else:
                                score += aligner.mismatch_score
                        else:
                            score += freq * aligner.substitution_matrix[char, prof_char]
        return score

    n_sequences = len(initial_alignment)
    if n_sequences <= 0:
        print('Error: iterative_refinement() received empty alignment')
        exit(1)
    # we need at least 3 sequences to do anything meaningful
    if n_sequences < 3:
        return initial_alignment
    current_alignment = initial_alignment
    current_score = sum(score_against_profile(seq, create_profile(current_alignment)) for seq in current_alignment)
    #
    for _ in range(max_iterations):
        new_alignment = []
        for i in range(len(current_alignment)):
            # remove the current sequence from the alignment
            temp_alignment = current_alignment[:i] + current_alignment[i+1:]
            # align the current sequence to a profile created from the remaining sequences
            profile = create_profile(temp_alignment)
            consensus = ''.join(max(col, key=col.get) for col in profile.values())
            alignment = aligner.align(consensus, current_alignment[i].replace('-', ''))[0]
            new_alignment.append(str(alignment[1]))
        max_len = max(len(seq) for seq in new_alignment)
        new_alignment = [seq.ljust(max_len, '-') for seq in new_alignment]
        #
        new_score = sum(score_against_profile(seq, create_profile(new_alignment)) for seq in new_alignment)
        if abs(current_score) < 1e-6:
            if new_score < current_score:
                break
        else:
            relative_improvement = (new_score - current_score) / abs(current_score)
            if relative_improvement < improvement_threshold:
                break
        current_alignment = new_alignment
        current_score = new_score
    #
    return remove_gap_columns(current_alignment)


def get_final_tvr_consensus(aligned, default_char=None, min_coverage=3, max_gap_frac=0.60, untrustworthy_chars='', tiebreak_adj={}):
    aligned_rstrip = [n.rstrip('-') for n in aligned]
    n_sequences = len(aligned_rstrip)
    if n_sequences <= 0:
        print('Error: get_final_tvr_consensus() received empty alignment')
        exit(1)
    # not enough sequences present, just pretend the first one is the consensus
    if n_sequences < min_coverage:
        return aligned_rstrip[0]
    consensus = []
    for col in zip_longest(*aligned_rstrip, fillvalue=' '):
        sorted_chars = Counter(col).most_common()
        char_count = defaultdict(int, {n[0]:n[1] for n in sorted_chars})
        # if not enough reads cover this position then report the default char
        coverage = n_sequences - char_count[' ']
        if coverage < min_coverage:
            if default_char is None:
                continue
            consensus.append(default_char)
            continue
        gap_count = char_count['-']
        sorted_chars = [n for n in sorted_chars if n[0] not in ' -']
        # if the most frequent character is noncanonical, trustworthy and has enough support, lets go with that
        if len(sorted_chars) >= 1 and sorted_chars[0][0] not in untrustworthy_chars and sorted_chars[0][1] >= min_coverage:
            consensus.append(sorted_chars[0][0])
            continue
        # ok, how about the second most frequent character?
        if len(sorted_chars) >= 2 and sorted_chars[1][0] not in untrustworthy_chars and sorted_chars[1][1] >= min_coverage:
            consensus.append(sorted_chars[1][0])
            continue
        # skip gaps
        if gap_count / n_sequences > max_gap_frac:
            continue
        # is it unanimous?
        if len(sorted_chars) == 1:
            consensus.append(sorted_chars[0][0])
            continue
        # are we untied?
        if sorted_chars[0][1] > sorted_chars[1][1]:
            consensus.append(sorted_chars[0][0])
            continue
        # ok we're tied for most frequent character and need to do a bit of cheating to pick the one we want
        adj_scores = []
        for (candidate_letter, _) in sorted_chars:
            if candidate_letter in tiebreak_adj:
                adj_scores.append((tiebreak_adj[candidate_letter], candidate_letter))
            else:
                adj_scores.append((0, candidate_letter))
        adj_scores = sorted(adj_scores, reverse=True)
        consensus.append(adj_scores[0][1])
    consensus = ''.join(consensus)
    return consensus


def get_nucl_consensus(sequences):
    my_dists = np.zeros((len(sequences), len(sequences)))
    aligner = get_aligner_object(gap_bool=(True,True), which_type='tvr')
    initial_msa = progressive_alignment(sequences, my_dists, aligner)
    refined_msa = iterative_refinement(initial_msa, aligner)
    consensus = get_final_tvr_consensus(refined_msa)
    return consensus


def quick_compare_tvrs(tvr_1, tvr_2):
    aligner = get_aligner_object(gap_bool=(True,False), which_type='tvr')
    return tvr_distance(tvr_1, tvr_2, aligner, adjust_lens=False, min_viable=False) / MAX_SEQ_DIST
