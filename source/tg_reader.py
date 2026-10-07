import pysam
import gzip
import sys
from collections import deque
from concurrent.futures import ProcessPoolExecutor
from contextlib import nullcontext

_BGZF_EOF = bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000')

#
# accepts fq / fq.gz / fa / fa.gz / bam / cram
#
# handles fa files with newlines in read sequence
#


class TG_Reader:
    def __init__(self, input_filename, replace_tabs_with_spaces=True, verbose=True,
                 ref_fasta='', bam_threads=1):
        self.replace_tabs_with_spaces = replace_tabs_with_spaces
        self.verbose = verbose
        fnl = input_filename.lower()
        #
        if fnl[-2:] == 'fq' or fnl[-5:] == 'fq.gz' or fnl[-5:] == 'fastq' or fnl[-8:] == 'fastq.gz':
            self.filetype = 'FASTQ'
        elif fnl[-2:] == 'fa' or fnl[-5:] == 'fa.gz' or fnl[-5:] == 'fasta' or fnl[-8:] == 'fasta.gz':
            self.filetype = 'FASTA'
        elif fnl[-3:] == 'bam':
            self.filetype = 'BAM'
        elif fnl[-4:] == 'cram':
            self.filetype = 'CRAM'
        else:
            print('Error: unknown file type given to TG_Reader():')
            print(' - acceptable input types: fq / fq.gz / fa / fa.gz / bam / cram')
            exit(1)
        #
        if self.filetype in ['BAM', 'CRAM']:
            if self.verbose:
                print('getting reads from ' + self.filetype + '...')
            if self.filetype == 'BAM':
                # pysam cannot combine threaded decoding with ignore_truncation.
                # Keep the original tolerant reader for BAMs without an EOF block.
                if bam_threads > 1:
                    try:
                        with open(input_filename, 'rb') as bam:
                            bam.seek(-len(_BGZF_EOF), 2)
                            if bam.read() != _BGZF_EOF:
                                bam_threads = 1
                    except OSError:
                        bam_threads = 1
                self.f = pysam.AlignmentFile(input_filename, "rb", ignore_truncation=bam_threads == 1,
                                             check_sq=False, threads=bam_threads)
            else:
                if ref_fasta == '': # cram without reference will almost certainly break, but try anyway
                    print()
                    print('Warning: trying to open a cram without a specified reference...')
                    print()
                    self.f = pysam.AlignmentFile(input_filename, "rc", ignore_truncation=True, check_sq=False)
                else:
                    self.f = pysam.AlignmentFile(input_filename, "rc", ignore_truncation=True, check_sq=False, reference_filename=ref_fasta)
            self.alns = self.f.fetch(until_eof=True)
        else:
            if fnl[-3:] == '.gz':
                if self.verbose:
                    print('getting reads from gzipped ' + self.filetype + '...')
                self.f = gzip.open(input_filename, 'rt')
            else:
                if self.verbose:
                    print('getting reads from ' + self.filetype + '...')
                self.f = open(input_filename, 'r')
        #
        self.buffer = []
        self.current_readname = None

    #
    # returns (readname, readsequence, qualitysequence, is_supplementary)
    #
    def get_next_read(self, with_quality=True):
        if self.filetype == 'FASTQ':
            my_name = self.f.readline().strip()[1:]
            if not my_name:
                return ('','','',False)
            if self.replace_tabs_with_spaces:
                my_name = my_name.replace('\t', ' ')
            my_read = self.f.readline().strip()
            _       = self.f.readline().strip()
            my_qual = self.f.readline().strip()
            if not with_quality:
                my_qual = ''
            return (my_name, my_read, my_qual, False)
        #
        elif self.filetype == 'FASTA':
            if self.current_readname is None:
                self.current_readname = self.f.readline().strip()[1:]
            if not self.current_readname:
                return ('','','',False)
            if self.replace_tabs_with_spaces:
                self.current_readname = self.current_readname.replace('\t', ' ')
            hit_eof = False
            while True:
                my_dat = self.f.readline().strip()
                if not my_dat:
                    hit_eof = True
                    break
                self.buffer.append(my_dat)
                if '>' in self.buffer[-1]:
                    break
            if hit_eof:
                out_dat = (self.current_readname, ''.join(self.buffer), '', False)
                self.current_readname = None
                self.buffer = []
            else:
                out_dat = (self.current_readname, ''.join(self.buffer[:-1]), '', False)
                self.current_readname = self.buffer[-1][1:]
                self.buffer = []
            return out_dat
        #
        elif self.filetype in ['BAM', 'CRAM']:
            try:
                aln = next(self.alns)
                return (aln.qname, aln.query_sequence,
                        aln.qual if with_quality else '', aln.is_supplementary)
            # we reached the end of file
            except StopIteration:
                return ('','','',False)
            # this can happen if file is truncated
            except OSError:
                return ('','','',False)

    #
    # returns list of [(readname1, readsequence1, qualitysequence1, issup1), (readname2, readsequence2, qualitysequence2, issup2), ...]
    #
    def get_all_reads(self):
        all_read_dat = []
        while True:
            read_dat = self.get_next_read()
            if not read_dat[0]:
                break
            all_read_dat.append((read_dat[0], read_dat[1], read_dat[2], read_dat[3]))
        return all_read_dat

    def close(self):
        self.f.close()


def quick_grab_all_reads(fn):
    #
    # a convenience function to grab all reads from a file in a single line
    #
    my_reader = TG_Reader(fn, verbose=False)
    all_read_dat = my_reader.get_all_reads()
    my_reader.close()
    return all_read_dat


def quick_grab_all_reads_nodup(fn, min_len=None):
    #
    # a modified version for ensuring no duplicates (e.g. reading in a bam with multimapped reads)
    #
    my_reader = TG_Reader(fn, verbose=False)
    all_read_dat = my_reader.get_all_reads()
    my_reader.close()
    by_readname = {}
    reads_filtered = 0
    for n in all_read_dat:
        if min_len is None or len(n[1]) >= min_len:
            by_readname[n[0]] = (n[1], n[2])
        else:
            reads_filtered += 1
    out_readdat = []
    for k in by_readname:
        out_readdat.append((k, by_readname[k][0], by_readname[k][1]))
    return (out_readdat, reads_filtered)


def _find_telomere_reads(sequences, kmer, reverse_kmer, min_hits):
    """Return one match flag per sequence, using the original non-overlapping count."""
    return [seq.count(kmer) >= min_hits or seq.count(reverse_kmer) >= min_hits
            for seq in sequences]


def _screen_and_compress(reads, kmer, reverse_kmer, min_hits):
    matches = _find_telomere_reads(
        [sequence for _, sequence in reads], kmer, reverse_kmer, min_hits)
    selected = ''.join(f'>{name}\n{sequence}\n' for (name, sequence), matched
                       in zip(reads, matches) if matched)
    # Concatenated gzip members can be read as one file by gzip and TG_Reader.
    compressed = gzip.compress(selected.encode(), compresslevel=6, mtime=0) if selected else b''
    return matches, compressed


def extract_telomere_reads(input_files, output_file, kmer, reverse_kmer,
                           min_hits, num_processes=1, ref_fasta=''):
    """Stream the initial repeat screen in bounded, ordered batches."""
    all_readcount = tel_readcount = sup_readcount = 0
    total_bp_all = total_bp_tel = 0
    readlens_all, readlens_tel = [], []
    batch = []
    batch_bp = 0
    pending = deque()
    max_pending = 2 * num_processes

    def record_batch(lengths, matches):
        nonlocal all_readcount, tel_readcount, total_bp_all, total_bp_tel
        for read_len, matched in zip(lengths, matches):
            all_readcount += 1
            total_bp_all += read_len
            readlens_all.append(read_len)
            if matched:
                tel_readcount += 1
                total_bp_tel += read_len
                readlens_tel.append(read_len)

    def submit_batch(executor, output):
        nonlocal batch, batch_bp
        if executor is None:
            matches = _find_telomere_reads(
                [read[1] for read in batch], kmer, reverse_kmer, min_hits)
            for (name, sequence), matched in zip(batch, matches):
                if matched:
                    output.write(f'>{name}\n{sequence}\n')
            record_batch([len(sequence) for _, sequence in batch], matches)
        else:
            pending.append(([len(sequence) for _, sequence in batch], executor.submit(
                _screen_and_compress, batch,
                kmer, reverse_kmer, min_hits)))
            if len(pending) >= max_pending:
                lengths, future = pending.popleft()
                matches, compressed = future.result()
                record_batch(lengths, matches)
                output.write(compressed)
        batch = []
        batch_bp = 0

    pool = ProcessPoolExecutor(max_workers=num_processes) if num_processes > 1 else nullcontext()
    output_file_handle = (open(output_file, 'wb') if num_processes > 1
                          else gzip.open(output_file, 'wt'))
    with pool as executor, output_file_handle as output:
        for input_file in input_files:
            reader = TG_Reader(input_file, verbose=False, ref_fasta=ref_fasta,
                               bam_threads=num_processes)
            try:
                while True:
                    name, sequence, _, is_supplementary = reader.get_next_read(with_quality=False)
                    if not name:
                        break
                    if not sequence:
                        continue
                    if is_supplementary:
                        sup_readcount += 1
                        continue
                    batch.append((name, sequence))
                    batch_bp += len(sequence)
                    if len(batch) >= 256 or batch_bp >= 4_000_000:
                        submit_batch(executor, output)
            finally:
                reader.close()
        if batch:
            submit_batch(executor, output)
        while pending:
            lengths, future = pending.popleft()
            matches, compressed = future.result()
            record_batch(lengths, matches)
            output.write(compressed)
        if executor is not None and tel_readcount == 0:
            output.write(gzip.compress(b'', mtime=0))

    return (all_readcount, tel_readcount, sup_readcount, total_bp_all,
            total_bp_tel, readlens_all, readlens_tel)


if __name__ == '__main__':
    #
    IN_READS_TEST = sys.argv[1]
    my_reader = TG_Reader(IN_READS_TEST)
    while True:
        read_dat = my_reader.get_next_read()
        if not read_dat[0]:
            break
        print(read_dat)
    my_reader.close()
