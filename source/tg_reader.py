import pysam
import gzip
import sys

#
# accepts fq / fq.gz / fa / fa.gz / bam / cram
#
# handles fa files with newlines in read sequence
#


class TG_Reader:
    def __init__(self, input_filename, replace_tabs_with_spaces=True, verbose=True, ref_fasta=''):
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
                self.f = pysam.AlignmentFile(input_filename, "rb", ignore_truncation=True, check_sq=False)
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
    def get_next_read(self):
        if self.filetype == 'FASTQ':
            my_name = self.f.readline().strip()[1:]
            if not my_name:
                return ('','','',False)
            if self.replace_tabs_with_spaces:
                my_name = my_name.replace('\t', ' ')
            my_read = self.f.readline().strip()
            _       = self.f.readline().strip()
            my_qual = self.f.readline().strip()
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
                return (aln.qname, aln.query_sequence, aln.qual, aln.is_supplementary)
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
