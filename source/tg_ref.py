import gzip
import os


class ReferenceFasta:
    def __init__(self, reference_fasta_fn):
        self.ref_fn   = reference_fasta_fn
        self.ref_inds = self.index_ref()
        self.ri_rev   = {self.ref_inds[n][0]:n for n in range(len(self.ref_inds))}

    def get_contig_names(self):
        return [n[0] for n in self.ref_inds]

    def get_refseq(self, contig_name):
        if contig_name not in self.ri_rev:
            print('Warning: contig not found in reference file:', contig_name)
            return None
        my_ri = self.ri_rev[contig_name]
        if self.ref_fn[-3:] == '.gz':
            # reading entire reference into memory, this is suboptimal but indexing gzipped files seems complicated
            with gzip.open(self.ref_fn, 'rt') as ref_file:
                fr = ref_file.read()
            my_dat = fr[int(self.ref_inds[my_ri][1]):int(self.ref_inds[my_ri][2])].replace('\n','')
        else:
            with open(self.ref_fn,'r') as ref_file:
                ref_file = open(self.ref_fn,'r')
                ref_file.seek(self.ref_inds[my_ri][1])
            my_dat = ref_file.read(int(self.ref_inds[my_ri][2]) - int(self.ref_inds[my_ri][1])).replace('\n','')
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
        if self.ref_fn[-3:] == '.gz':
            fn = None
        #
        ref_inds = []
        if fn is not None:
            with open(fn,'r') as fai:
                for line in fai:
                    splt = line[:-1].split('\t')
                    seqlen  = int(splt[1])
                    offset  = int(splt[2])
                    linelen = int(splt[3])
                    nlines  = seqlen // linelen
                    if seqlen % linelen != 0:
                        nlines += 1
                    ref_inds.append((splt[0], offset, offset + seqlen + nlines, seqlen))
        else:
            if verbose:
                print('fasta index not found, creating one...')
            prevR = None
            prevP = None
            seqlen = 0
            if self.ref_fn[-3:] == '.gz':
                f = gzip.open(self.ref_fn, 'rt')
            else:
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
