# Liz_Sanger_seq

This respository contains rawdata, R routines and results of the analysis of Sanger sequencing data from different samples.

## Sequencing data

Raw sequencing data is based on Sanger sequencing targeting exons of gene *VWF* (Von Willebrand factor). Specifically, we have

  * sequencing data of proband from 2016 for all 52 exons of gene *VWF*, and
  * sequencing data of father, mother, brother and sister from 2018 for exons 23, 27, 33, 39 of gene *VWF*.

Filenames have the following format:

````
<source>-<id1>-<exon_no><primer_orientation>_<id2>.ab1
````

Example 1:

`1-VMF-23R_E06.ab1`:

  * `<source> = 1`
  * `<id1> = VMF`
  * `<exon_no> = 23`
  * `<primer_orientation> = R`
  * `<id2> = E06`

Example 2:

`QJL-VWF-50-VWF-50F_C07.ab1`

  * `<source> = QJL`
  * `<id1> = VMF-50-VWF`
  * `<exon_no> = 50`
  * `<primer_orientation> = F`
  * `<id2> = C07`

Note that primer orientation can change for the same exon across different sources. For example we have the following primer orientations for exons 23, 27, 33, 39:


Exons | QJL (Proband) | 1 (Father) | 2 (Mother) | 3 (Sister) | 4 (Brother)
------ | ------- | ------ | ------ | ------ | -------
E23 | Reverse | Reverse | Reverse | Reverse | Reverse
E27 | Forward | Forward | Forward | Forward | Forward
E33 | Forward | Forward | Forward | Reverse | Forward
E39 | Forward | Forward | Forward | Forward | Forward


## *VWF* Reference sequence

We use [*NM_000552.4*](https://www.ncbi.nlm.nih.gov/nuccore/NM_000552) from NCBI's RefSeq database as reference sequence.

Note that [*NM_000552.3*](https://www.ncbi.nlm.nih.gov/nuccore/NM_000552.3) is missing the first 5 nucleotides but is otherwise identical to *NM_000552.4*. Results from a pairwise alignment analysis are given in file [aln-NM_000052.3-NM_000552.4.txt](aln-NM_000052.3-NM_000552.4.txt).


## Multiple sequence alignment

We perform multiple sequence alignment of the reference sequence and primary call Sanger sequences for proband, father, mother, sister and brother, and for every exon using the R library [`msa`](https://bioconductor.org/packages/release/bioc/html/msa.html).

  1. Primary/secondary basecall traceplots for all samples are given in `./plots/<sample>.pdf`, where `<sample>` is the filename as detailed above without the `.ab1` extension.  
  2. MSA results are given in `./plots/<exon_no>.pdf`, where `<exon_no>` is the exon number. Note that father, mother, sister, brother data is only available for exons 23, 27, 33, 39.


## R library requirements

R script `script.R` was tested using the following libraries:

  * `Biostrings_2.46.0`
  * `GenomicRanges_1.30.3`
  * `msa_1.10.0`
  * `sangerSeqR_1.14.0`
  * `tidyverse_1.2.1`
