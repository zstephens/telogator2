import sys
import subprocess

from source.tg_reader import TG_Reader
from source.tg_util import rm

MUSCLE_VERS_STRING = 'MUSCLE v3.8'
MUSCLE_URL         = 'https://drive5.com/muscle/downloads_v3.htm'

#
#
#
def check_muscle_version(MUSCLE_EXE):
    sys.stdout.write('checking MUSCLE version...')
    sys.stdout.flush()
    found_muscle = False
    try:
        output = subprocess.check_output([MUSCLE_EXE, '-version']).decode('utf-8')
        if output[:len(MUSCLE_VERS_STRING)] == MUSCLE_VERS_STRING:
            found_muscle = True
    except FileNotFoundError:
        pass
    except PermissionError:
        pass
    if found_muscle:
        sys.stdout.write(' ok\n')
        sys.stdout.flush()
    else:
        sys.stdout.write(' failed\n')
        sys.stdout.flush()
        print('Error: MUSCLE v3.8 exe not found')
        print(' - MUSCLE can be downloaded at:', MUSCLE_URL)
        exit(1)

#
#
#
def get_muscle_msa(input_sequences, muscle_exe, working_dir='', char_score_adj={}, max_gap_frac=0.60, noncanon_cheat=None, mode='amino'):
    # check if we received all empty strings for some reason
    if len([n for n in input_sequences if len(n) > 0]) == 0:
        return ['', '']
    # write sequences to a temp fasta
    temp_fasta = working_dir + 'clust_sequences.fa'
    f = open(temp_fasta, 'w')
    for i in range(len(input_sequences)):
        if len(input_sequences[i]) > 0:
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
    score_param += ' -maxiters 2'   # use if muscle is crashing on "refining bipartite" steps
    cmd = muscle_exe + ' -in ' + temp_fasta + ' -out ' + aln_fasta + ' ' + score_param
    #
    try:
        output = subprocess.check_output(cmd.split(' '), stderr=subprocess.STDOUT, text=True)
    except subprocess.CalledProcessError as exc:
        print('Error: MUSCLE returned an error:', exc.returncode)
        print(exc.output)
        exit(1)
    #
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
    ####    print(n)
    ####exit(1)
    # cleanup
    rm(temp_fasta)
    rm(aln_fasta)
    rm(muscle_log)
    rm(matrix)
    # get consensus
    consensus_seq = []
    for i in range(len(out_seq[0])):    # we're gonna hope that muscle worked as expected and all seq are same len
        char_count = {}
        gap_count  = 0
        for j in range(len(out_seq)):
            if out_seq[j][i] == '-':
                gap_count += 1
            else:
                if out_seq[j][i] not in char_count:
                    char_count[out_seq[j][i]] = 0
                char_count[out_seq[j][i]] += 1
        # if first or second most frequent character is non-canonical (and has enough support), lets just go with that
        if noncanon_cheat is not None:
            (canon_char, min_noncanon_reads) = noncanon_cheat
            sorted_chars = sorted([(char_count[k],k) for k in char_count.keys()], reverse=True)
            if sorted_chars[0][0] >= min_noncanon_reads and sorted_chars[0][1] != canon_char:
                consensus_seq.append(sorted_chars[0][1])
                continue
            if len(sorted_chars) >= 2 and sorted_chars[1][0] >= min_noncanon_reads and sorted_chars[1][1] != canon_char:
                consensus_seq.append(sorted_chars[1][1])
                continue
        #
        if float(gap_count)/len(out_seq) > max_gap_frac:
            continue
        #
        candidates = [(char_count[k],k) for k in char_count.keys() if char_count[k] == max(char_count.values())]
        if len(candidates) == 1:
            consensus_seq.append(candidates[0][1])
        else:   # tie-breaking logic
            adj_scores = []
            for candidate in candidates:
                if candidate[1] in char_score_adj:
                    adj_scores.append((char_score_adj[candidate[1]], candidate[1]))
                else:
                    adj_scores.append((0, candidate[1]))
            adj_scores = sorted(adj_scores, reverse=True)
            consensus_seq.append(adj_scores[0][1])
    consensus_seq = ''.join(consensus_seq)
    #
    return [out_seq, consensus_seq]

#
#
#
def write_amino_scoring_matrix(fn):
    f = open(fn, 'w')
    f.write('   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *' + '\n')
    f.write('A  2  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -4' + '\n')
    f.write('R  0  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
    f.write('N  0 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
    f.write('D  0 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
    f.write('C  1 -4 -4 -4  0 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
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

#
#
#
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
