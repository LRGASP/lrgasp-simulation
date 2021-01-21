# A pipeline for simulating RNA sequencing data


## About the pipeline

This pipeline is designed for simulating RNA sequencing datasets that includes
PacBio, Oxford Nanopore and Illumina reads based on the provided reference data and 
expression profile. Basically, it is a wrapper for several simulation tools: 
- [IsoSeqSim](https://github.com/yunhaowang/IsoSeqSim)
- [Trans-NanoSim](https://github.com/bcgsc/NanoSim)
- [RSEM simulator](http://deweylab.biostat.wisc.edu/rsem/README.html)


## Getting expression profile

To create an expression profile, you need to estimate transcript abundances 
using real long-read sequencing data. To do so, add [`minimap2`](https://github.com/lh3/minimap2) to your
`$PATH` variable and run

``` quantify.py -r <REFERENCE_TRANSCRIPTS.fasta> --fastq <READS.fastq -o COUNTS.tsv```

Available options are:

``` --reference_transcripts, -r``` reference transcriptome in FASTA format

``` --fastq, -f``` long RNA reads in FASTQ format

``` --output, -o``` output file with abundances (counts and TPM) in TSV format

``` --threads, -t``` number of threads for `minimap2`

``` --mandatory, -m``` file with a list of mandatory transcripts to be included,
                       counts are assigned randomly, can be a TSV with transcript ids in the first column;
                       this option is used to provide artificial "novel" transcripts;
                       make sure they are included in the reference

``` --seed, -m``` randomizer seed

## Simulating reads

To simulate reads run

``` simulate.py --reference_dir <PATH/TO/REFERENCES/> --reference_name <REFERENCE_NAME> --counts <COUNTS.tsv> -o <OUTPUT_DIR> ```

For example, to run on test data included in the repository launch

``` simulate.py --reference_dir data/ --reference_name test --counts data/test.counts.tsv -o test_simulation ```

Available options are:

``` --reference_dir, -r``` folder with reference data (must contain genome, transcriptome and annotation)

``` --reference_name, -n``` prefix of reference files (files are `.genome.fasta`, `.transcripts.fasta` and `.annotation.gtf`) 

``` --output, -o``` output folder

``` --threads, -t``` number of threads