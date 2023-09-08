import gzip
import os


class ReferenceFasta:
    def __init__(self, reference_fasta_fn):
        self.ref_fn   = reference_fasta_fn
        self.ref_inds = self.index_ref()
        self.ri_rev   = {self.ref_inds[n][0]:n for n in range(len(self.ref_inds))}

    def get_refseq(self, contig_name):
        if contig_name not in self.ri_rev:
            print('Warning: contig not found in reference file:', contig_name)
            return None
        my_ri = self.ri_rev[contig_name]
        ref_file = open(self.ref_fn,'r')
        ref_file.seek(self.ref_inds[my_ri][1])
        my_dat = ''.join(ref_file.read(int(self.ref_inds[my_ri][2]) - int(self.ref_inds[my_ri][1])).split('\n'))
        ref_file.close()
        return my_dat

    def index_ref(self, verbose=False):
        fn = None
        if os.path.isfile(self.ref_fn + 'i'):
            if verbose:
                print('found fasta index ' + self.ref_fn + 'i')
            fn = self.ref_fn + 'i'
        if os.path.isfile(self.ref_fn + '.fai'):
            if verbose:
                print('found fasta index ' + self.ref_fn + '.fai')
            fn = self.ref_fn + '.fai'
        #
        ref_inds = []
        if fn is not None:
            fai = open(fn,'r')
            for line in fai:
                splt = line[:-1].split('\t')
                seqlen  = int(splt[1])
                offset  = int(splt[2])
                linelen = int(splt[3])
                nlines  = seqlen // linelen
                if seqlen % linelen != 0:
                    nlines += 1
                ref_inds.append((splt[0], offset, offset + seqlen + nlines, seqlen))
            fai.close()
        else:
            if verbose:
                print('fasta index not found, creating one...')
            prevR = None
            prevP = None
            seqlen = 0
            f = open(self.ref_fn, 'r')
            while True:
                data = f.readline()
                if not data:
                    ref_inds.append((prevR, prevP, f.tell() - len(data), seqlen))
                    break
                if data[0] == '>':
                    if prevP is not None:
                        ref_inds.append((prevR, prevP, f.tell() - len(data), seqlen))
                    seqlen = 0
                    prevP  = f.tell()
                    prevR  = data[1:-1]
                else:
                    seqlen += len(data) - 1
            f.close()
        return ref_inds
