# Telogator2
A method for measuring allele-specific TL and characterizing telomere variant repeat (TVR) sequences from long reads. Telogator2 was designed to work with long reads that are capable of representing telomere regions without high sequencing error rates, and has been tested primarily on PacBio HiFi reads and ONT reads basecalled with Dorado.

If this software has been useful for your work, please cite us at:

(Publication currently in submission)

In the mean time, see our preprint:

https://www.researchsquare.com/article/rs-3673861/v1


## Dependencies:

Telogator2 dependencies can be easily installed via conda:

```bash
# create conda environment
conda env create -f conda_env_telogator2.yaml

# activate environment
conda activate telogator2
```


## Running Telogator2:

Telogator2 can be run in single command:

```bash
python telogator2.py -i input.fq \ 
                     -o results/ \ 
                     --muscle /path/to/muscle \ 
                     --minimap2 /path/to/minimap2
```

`-i` accepts fa / fa.gz / fq / fq.gz / bam. Multiple arguments can be provided (e.g. `-i reads1.fa reads2.fa`).

Telogator2 requires that [muscle v3.8](https://drive5.com/muscle/downloads_v3.htm) is installed in order to produce consensus TVR sequences. If the dependencies were installed via conda then muscle should be found in the `envs/telogator2/` directory. Additionally, Telogator2 requires a path to an aligner, via either the `--minimap2`, `--winnowmap`, or `--pbmm2` input parameters.



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



## Test data

Telomere reads for HG002 can be found in the `test_data/` directory.

```
HiFi reads (~70x): hg002-telreads_pacbio.fa.gz
ONT reads  (~25x): hg002-telreads_ont.fa.gz
```



## Telogator reference

The reference sequence used for telomere anchoring currently contains the first and last 500kb of sequences from the following T2T assemblies:

* `T2T-chm13` - https://github.com/marbl/CHM13
* `T2T-yao` - https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA017932
* `T2T-cn1` - https://github.com/T2T-CN1/CN1
* `T2T-hg002` - https://github.com/marbl/hg002
* `T2T-ksa001` - https://github.com/bio-ontology-research-group/KSA001
* `T2T-i002c` - https://github.com/LHG-GG/I002C

More will be added as they become available.
