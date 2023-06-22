import pysam
import gzip
import sys

#
# accepts fq / fq.gz / fa / fa.gz / bam
#
# handles fa files with newlines in read sequence
#


class TG_Reader:
    def __init__(self, input_filename, replace_tabs_with_spaces=True, verbose=True):
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
        else:
            print('Error: unknown file suffix given to TG_Reader:')
            print(input_filename)
            exit(1)
        #
        if self.filetype == 'BAM':
            if self.verbose:
                print('getting reads from ' + self.filetype + '...')
            self.f = pysam.AlignmentFile(input_filename, "rb", ignore_truncation=True, check_sq=False)
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
    # returns (readname, readsequence, qualitysequence)
    #
    def get_next_read(self):
        if self.filetype == 'FASTQ':
            my_name = self.f.readline().strip()[1:]
            if not my_name:
                return ('','','')
            if self.replace_tabs_with_spaces:
                my_name = my_name.replace('\t', ' ')
            my_read = self.f.readline().strip()
            skip    = self.f.readline().strip()
            my_qual = self.f.readline().strip()
            return (my_name, my_read, my_qual)
        #
        elif self.filetype == 'FASTA':
            if self.current_readname is None:
                self.current_readname = self.f.readline().strip()[1:]
            if not self.current_readname:
                return ('','','')
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
                out_dat = (self.current_readname, ''.join(self.buffer), '')
                self.current_readname = None
                self.buffer = []
            else:
                out_dat = (self.current_readname, ''.join(self.buffer[:-1]), '')
                self.current_readname = self.buffer[-1][1:]
                self.buffer = []
            return out_dat
        #
        elif self.filetype == 'BAM':
            try:
                aln = next(self.alns)
                my_name = aln.qname
                try:
                    my_np = aln.get_tag('np')
                    my_rq = aln.get_tag('rq')
                    my_name += ' np={0} rq={1:0.6f}'.format(my_np, my_rq)
                # np or rq tag is not present
                except KeyError:
                    pass
                # get read sequence directly from SAM entry instead of using aln.query (which doesn't include softclipped bases)
                aln_readdat = str(aln).split('\t')[9]
                return (my_name, aln_readdat, aln.qual)
            # this can happen if file is truncated
            except OSError:
                return ('','','')

    #
    # returns list of [(readname1, readsequence1, qualitysequence1), (readname2, readsequence2, qualitysequence2), ...]
    #
    def get_all_reads(self):
        all_read_dat = []
        while True:
            read_dat = self.get_next_read()
            if not read_dat[0]:
                break
            all_read_dat.append((read_dat[0], read_dat[1], read_dat[2]))
        return all_read_dat

    def close(self):
        self.f.close()

#
# a convenience function to grab all reads from a file in a single line
#
def quick_grab_all_reads(fn):
    my_reader = TG_Reader(fn, verbose=False)
    all_read_dat = my_reader.get_all_reads()
    my_reader.close()
    return all_read_dat

#
# a quick test
#


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
