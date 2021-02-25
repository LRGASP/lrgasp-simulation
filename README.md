# A pipeline for simulating RNA sequencing data


## About the pipeline

This pipeline is designed for simulating RNA sequencing datasets that include
PacBio, Oxford Nanopore and Illumina reads based on the provided reference data and 
expression profile. Basically, it is a wrapper for several simulation tools: 
- [IsoSeqSim](https://github.com/yunhaowang/IsoSeqSim)
- [Trans-NanoSim](https://github.com/bcgsc/NanoSim)
- [RSEM simulator](http://deweylab.biostat.wisc.edu/rsem/README.html)

The pipeline consists of main 3 steps:
- Preparing reference data, which includes insertting artificial novel isoforms.
  These "novel" isoforms can be obtained by mapping any mammalian reference transcripts 
  onto your genome (human or mouse) and processing them with [SQANTI3](https://github.com/ConesaLab/SQANTI3). 
- Quantifying transcript abundance, which requires any real long-read RNA dataset.
- Generating simulated reads and supplementary information, which will be available to the evaluators only.


### Requirements

- biopython
- gffutils
- pysam
- minimap2

## Preparing reference data

To prepare reference data for simulation you will need to obtain reference genome and
reference transcripts in FASTA format, and gene annotation in GFF/GTF format.
Furthermore, you will need to run [SQANTI3](https://github.com/ConesaLab/SQANTI3) on
reference transcripts from any organism of your choice (e.g. rat).
The isoforms can be then selected randomly or manually. If you want to
simulate reads based on reference transcripts only, simply omit `--sqanti_prefix` option.

To prepare reference data run:

``` prepare_reference_data.py -a <annotation.gtf> -t <transcripts.fa> -g <genome.fa> -q <SQANTI output prefix> --n_random_isoforms <int> -o <output_prefix>  ```

Available options are:

``` --output, -o ``` output prefix

```--reference_annotation, -a``` reference annotation (GTF/.db)

```--reference_transcripts, -t``` reference transcripts in FASTA format

```--reference_genome, -g``` reference genome in FASTA format

```--sqanti_prefix, -q``` prefix of SQANTI output (`_classification.txt` and `_corrected.gtf` are needed)

```--n_random_isoforms, -n ``` insert this number of random novel artificial isoforms into the annotation

```--isoform_list, -l``` insert only novel artificial isoforms from a given file

```--seed, -s``` randomizer seed [11]

If both `--n_random_isoforms` and `--isoform_list` are ignored, all isoforms reported by SQANTI will be inserted.

If you want to skip inserting novel artificial isoforms and simulate reads based on reference transcripts only, 
simply omit `--sqanti_prefix` option.

## Quantifying transcript abundance

To create an expression profile, you need to estimate transcript abundances 
using real long-read sequencing data. To do so, add [`minimap2`](https://github.com/lh3/minimap2) to your
`$PATH` variable and run

``` quantify.py -r <trnascripts.fasta> --fastq <reads.fastq> -o counts.tsv```

We recommend to use transcript sequences obtained at the previous step.

Available options are:

``` --reference_transcripts, -r``` reference transcriptome in FASTA format

``` --fastq, -f``` long RNA reads in FASTQ format

``` --output, -o``` output file with abundances (counts and TPM) in TSV format

``` --threads, -t``` number of threads for `minimap2`

``` --mandatory, -m``` file with a list of mandatory transcripts to be included,
                       counts are assigned randomly, can be a TSV with transcript ids in the first column;
                       this option is used to provide artificial "novel" transcripts;
                       make sure they are included in the reference;
                       we recommend to use a list of novel isoforms generated at the previous step;

``` --seed, -m``` randomizer seed

## Simulating reads

To simulate reads run

``` simulate.py --reference_dir <path/to/references/> --reference_name <reference_prefix> --counts <counts.tsv> -o <output_dir> ```

For example, to run on test data included in the repository launch

``` simulate.py --reference_dir data/ --reference_name test --counts data/test.counts.tsv -o test_simulation ```

Available options are:

``` --reference_prefix, -r``` prefix for reference files (files are `.genome.fasta`, `.transcripts.fasta` and `.annotation.gtf`);
                              use the output of `prepare_reference_data.py`;

```--counts, -c``` transcript abundances in TSV format (output of `quantify.py`)

``` --output, -o``` output folder

``` --threads, -t``` number of threads


## Example 

Example data (Human chromosome 22) can be downloaded from:

section in progress

## Output

## Reference data

Human ([gencode v37](https://www.gencodegenes.org/human/)):
- ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz
- ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.transcripts.fa.gz

Mouse ([gencode M26](https://www.gencodegenes.org/mouse/)):
- ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/gencode.vM26.annotation.gtf.gz
- ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/gencode.vM26.transcripts.fa.gz