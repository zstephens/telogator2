# Telogator2
A method for measuring allele-specific TL and characterizing telomere variant repeat (TVR) sequences from long reads. Telogator2 was designed to work with long reads that are capable of representing telomere regions without high sequencing error rates, and has been tested primarily on PacBio HiFi reads.

If this software has been useful for your work, please cite us at:

(Publication currently in submission)

Alternately, see our paper for the previous version of Telogator:

Stephens, Zachary, et al. "Telogator: a method for reporting chromosome-specific telomere lengths from long reads." *Bioinformatics* 38.7 (2022): 1788-1793. https://doi.org/10.1093/bioinformatics/btac005


## Dependencies:

Telogator2 dependencies can be easily installed via conda:

```bash
# create conda environment
conda env create -f conda_env_telogator2.yaml

# activate environment
conda activate telogator2
```



## (1) making T2T + alt subtel reference:

Obtain t2t reference sequence from https://github.com/marbl/CHM13 (`chm13v2.0.fa`)

Append alternate subtelomere assemblies:

```bash
cat chm13v2.0.fa resources/stong_subtels.fa > t2t-and-subtel.fa
samtools faidx t2t-and-subtel.fa
```


## (2) align reads to whole genome reference:

For PacBio CLR reads:

```bash
pbmm2 align t2t-and-subtel.fa clr-reads.fa.gz aln.bam --preset SUBREAD --sort
```

For PacBio HiFi reads:

```bash
pbmm2 align t2t-and-subtel.fa hifi-reads.fa.gz aln.bam --preset HiFi --sort
```

For Oxford Nanopore reads:

```bash
minimap2 -ax map-ont -t 6 -N 5 -Y -L -o aln-unsort.sam t2t-and-subtel.fa ont-reads.fa.gz
samtools sort -o aln.bam aln-unsort.sam
samtools index aln.bam
```


## (3) extract subset of reads that were mapped to subtelomere regions:

```bash
python3 get_subtel_reads.py -b aln.bam -i clr-reads.fa.gz -o subtel-reads.fa.gz
```


## (4) align read subset to telogator reference:

We recommend using the [winnowmap](https://github.com/marbl/Winnowmap) aligner for best results. An example for HiFi reads:

```bash
winnowmap -ax map-hifi -Y -W resources/repetitive-k15.txt -o subtel_aln-unsort.sam t2t-telogator-ref.fa subtel-reads.fa.gz
samtools sort -o subtel_aln.bam subtel_aln-unsort.sam
samtools index subtel_aln.bam
```

Though this alignment could be done using the same tools as in step 2, if desired.



## (5) run telogator2 on subtel-only alignment:

```bash
python3 telogator.py -i subtel_aln.bam -o telogator_out/
```

The `-p` input option specifies the number of processes to use (default: 4).

In order to save time during reprocessing, the `-i` input option can also accept `tel-data.p` files so that BAM files don't need to be re-parsed each time.

Telogator2 requires that the muscle multiple sequence aligner v3.8 is installed in order to produce consensus TVR sequences: https://drive5.com/muscle/downloads_v3.htm The path to the muscle executable is specified via the `-m` input option.



## Test data

To quickly test the functionality of Telogator2, we provided a subset of hg002 telomere reads in test_data/

```bash
python3 telogator.py -i test_data/hg002_chr1q.p -o telogator_out/
```

With the default number of processes (4) this command should complete in about 1 minute.



## Output files

The primary output files are:

* `tlens_by_chr.tsv` chromosome-specific telomere lengths
* `tlens_by_allele.tsv` allele-specific telomere lengths
* `tvr_clustering/tvr-consensus-*` plots of consensus sequence for each allele at each anchor position
* `tvr_clustering/tvr-reads-*` pileup of supporting reads at each anchor position

The main results are in `tlens_by_allele.tsv`, which has the following columns:

* `chr` anchor chromosome arm
* `position` anchor coordinate
* `allele_num` an ID number used for keeping track of alleles internally
* `allele_id` ID number for this specific TVR. Multimapped alleles will have the same value
* `TL_max` allele-specific telomere length
* `read_TLs` TL of each supporting read in the cluster
* `read_lengths` length of each read in the cluster
* `read_mapq` mapping quality of each read in the cluster
* `tvr_len` length of the cluster's TVR region
* `tvr_consensus` consensus TVR region sequence
* `supporting_reads` readnames of each read in the cluster
* `original_chr` chromosome arm that reads were originally mapped to prior to being grouped together (only applicable to clusters with 'blank' TVR regions)
