# Telogator2
A method for measuring allele-specific TL and characterizing telomere variant repeat (TVR) sequences from long reads.

If this software has been useful for your work, please cite us at:

Stephens, Z., & Kocher, J. P. (2024). Characterization of telomere variant repeats using long reads enables allele-specific telomere length estimation. *BMC bioinformatics*, 25(1), 194.

https://link.springer.com/article/10.1186/s12859-024-05807-5



## Dependencies:

Telogator2 dependencies can be easily installed via conda:

```bash
# create conda environment
conda env create -f conda_env_telogator2.yaml

# activate environment
conda activate telogator2
```


## Running Telogator2:

```bash
python telogator2.py -i input.fq \ 
                     -o results/ \ 
                     --minimap2 /path/to/minimap2
```

`-i` accepts fa, fa.gz, fq, fq.gz, or bam (multiple can be provided, e.g. `-i reads1.fa reads2.fa`). For Revio reads sequenced with SMRTLink13 and onward, we advise including both the "hifi" BAM and "fail" BAM as input to Telogator2.

An aligner executable must be specified, via either `--minimap2`, `--winnowmap`, or `--pbmm2`.



## Recommended settings

Sequencing platforms have different sequencing error types, as such we recommend running Telogator2 with different options based on which platform was used:

**PacBio Revio HiFi (30x)** - `-r hifi -n 4`  
**PacBio Sequel II (10x)** - `-r hifi -n 3`  
**Nanopore R10 (30x)** - `-r ont -n 4`  

Telogator2 may be unable to analyze older Nanopore data, as reads basecalled with Guppy have prohibitively high sequencing error rates in telomere regions.

For large datasets, such as data from enrichment methods described by [Karimian et al.](https://www.science.org/doi/abs/10.1126/science.ado0431) or [Schmidt et al.](https://www.nature.com/articles/s41467-024-48917-7), higher thresholds may be needed to reduce false positives: `-r ont -n 10`.  

By default Telogator2 is run with 4 processes. Runtime can be greatly reduced by specifying more, e.g. `-p 8` or `-p 16`, based on your system's available CPU resources.



## Test data

Telomere reads for HG002 can be found in the `test_data/` directory. These are full-sized datasets and may take several hours to run.

```
HiFi reads (~70x): hg002-telreads_pacbio.fa.gz
ONT reads  (~25x): hg002-telreads_ont.fa.gz
```

A smaller dataset is also provided, which should take no more than a couple minutes to complete:

```bash
python telogator2.py -i test_data/hg002-ont-1p.fa.gz \ 
                     -o results/ \ 
                     -r ont
```


## Output files

The primary output files are:

* `tlens_by_allele.tsv` allele-specific telomere lengths
* `all_final_alleles.png` plots of all alleles (TVR + telomere regions)
* `violin_atl.png` violin plot of ATLs at each chromosome arm

The main results are in `tlens_by_allele.tsv`, which has the following columns:

* `chr` anchor chromosome arm
  * subtelomeres that could not be aligned are labeled `chrU` for 'unmapped'
* `position` anchor coordinate
* `ref_samp` the specific T2T reference contig to which the subtelomere was aligned
* `allele_id` ID number for this specific allele
  * ids ending in `i` indicate subtelomeres that were aligned to known interstitial telomere regions. These alleles should likely be excluded from subsequent analyses.
* `TL_p75` ATL (reports 75th percentile by default)
* `read_TLs` ATL of each supporting read in the cluster
* `read_lengths` length of each read in the cluster
* `read_mapq` mapping quality of each read in the cluster
* `tvr_len` length of the cluster's TVR region
* `tvr_consensus` consensus TVR region sequence
* `supporting_reads` readnames of each read in the cluster



## Telogator reference

The reference sequence used for telomere anchoring currently contains the first and last 500kb of sequences from the following T2T assemblies:

* `T2T-chm13` - https://github.com/marbl/CHM13
* `T2T-yao` - https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA017932
* `T2T-cn1` - https://github.com/T2T-CN1/CN1
* `T2T-hg002` - https://github.com/marbl/hg002
* `T2T-ksa001` - https://github.com/bio-ontology-research-group/KSA001
* `T2T-i002c` - https://github.com/LHG-GG/I002C

More will be added as they become available.

## Mouse reference

Experimental support has been added for the [T2T-mouse genome](https://github.com/yulab-ql/mhaESC_genome):

```bash
python telogator2.py -i input.fq \ 
                     -o results/ \ 
                     -t resources/telogator-ref-mouse.fa.gz \ 
                     --minimap2 /path/to/minimap2
```

Note that if you choose winnowmap as the aligner that you will also need to add `--winnowmap-k15 resources/telogator-ref-mouse-k15.txt`.