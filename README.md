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


## (1) extract reads with telomere repeats:


```bash
python3 get_reads_with_kmer_hits.py -i input.fq.gz -o telomere-reads.fa.gz
```

With default parameters this script extracts reads that contain at least 5 TTAGGTTAGG (or CCTAACCTAA) repeats. `-i` accepts .fa / .fa.gz / .fq / .fq.gz / .bam. Multiple input files can be provided. `-o` can be .fa / .fa.gz / .fq / .fq.gz.


## (2) align read subset to telogator reference:

We recommend using the [winnowmap](https://github.com/marbl/Winnowmap) aligner for best results. An example for HiFi reads:

```bash
winnowmap -ax map-pb -Y -W resources/repetitive-k15.txt -o aln-unsort.sam resources/t2t-telogator-ref.fa telomere-reads.fa.gz
samtools sort -o aln.bam aln-unsort.sam
samtools index aln.bam
```

Though you could also use pbmm2, minimap2, or another aligner of your choice.


## (3) run telogator2:

```bash
python3 telogator2.py -i aln.bam -o telogator_out/
```

The `-p` input option specifies the number of processes to use (default: 4).

Telogator2 requires that the muscle multiple sequence aligner v3.8 is installed in order to produce consensus TVR sequences: https://drive5.com/muscle/downloads_v3.htm

The path to the muscle executable is specified via the `-m` input option. If the dependencies were installed via conda then muscle should be found in the `envs/telogator2/` directory.


## Test data

To quickly test the functionality of Telogator2, we provided a subset of hg002 telomere reads in test_data/

```bash
python3 telogator2.py -i test_data/hg002_chr1q.p -o telogator_out/
```


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



## Temporary files

During TVR clustering, intermediary files are produced in the `tvr_clustering/temp/` directory which may be useful for debugging purposes. E.g. if reads at a particular chromosome arm are not being clustered together in the expected manner, consider investigating the `dendrogram-*` files to see if the default value for cutting clusters should be adjusted. This parameter can be adjusted via the `-att` input option.
