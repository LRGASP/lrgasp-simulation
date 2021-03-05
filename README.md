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

- python3
- biopython
- gffutils
- pysam
- htseq
- numpy
- scipy
- pandas
- [NanoSim requirements](https://github.com/bcgsc/NanoSim)
- [RSEM requirements](https://github.com/deweylab/RSEM#compilation) (includes C++, R, perl)  
- minimap2

## Preparing reference data

To prepare reference data for simulation you will need to obtain reference genome and
reference transcripts in FASTA format, and gene annotation in GFF/GTF format.
Furthermore, you will need to run [SQANTI3](https://github.com/ConesaLab/SQANTI3) on
reference transcripts from any organism of your choice (e.g. rat).
The isoforms can be then selected randomly or manually. If you want to
simulate reads based on reference transcripts only, simply omit `--sqanti_prefix` option.

To prepare reference data run:

```bash
prepare_reference_data.py \
  -a <annotation.gtf> \
  -t <transcripts.fa> \
  -g <genome.fa> \
  -q <SQANTI output prefix> \
  --n_random_isoforms <int> \
  -o <output_prefix>
```

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

```bash
simulate.py \
  --reference_prefix  <path/to/references/prefix> \
  -o <output_dir>
```

For example, to run on test data included in the repository launch:

```bash
python simulate.py \
  --reference_prefix data/test_data/test \
  --counts data/test_data/test.counts.tsv \
  --test_mode \
  -o test_simulation
```

or simply run the premade script:
```bash
bash run_simulate_test.sh
```

Available options are:

``` --reference_prefix, -r``` prefix for reference files (files are `.genome.fasta`, `.transcripts.fasta` and `.annotation.gtf`);
                              use the output of `prepare_reference_data.py`;

```--counts, -c``` transcript abundances in TSV format (output of `quantify.py`)

``` --output, -o``` output folder

``` --threads, -t``` number of threads

``` --seed, -s``` random seed to use

``` --ont_type``` type of molecule to simulate, 'dRNA' or 'cDNA'

``` --illumina_count``` number of Illumina read pairs to simulate

``` --pb_count``` number of PacBio reads to simulate

``` --ont_count``` number of ONT reads to simulate


## Example

Example data (Human chromosome 22) can be found in `data/human_chr22.tar.gz`:

Unpack the data by running `tar -xzf human_chr22.tar.gz` in `data` folder and launch the following commands.

Step 1: prepare reference data with 50 artificial novel isoforms

```prepare_reference_data.py --reference_annotation data/human_chr22/gencode.v32.annotation.chr22.gtf --reference_transcripts data/human_chr22/gencode.v32.transcripts.chr22.fa --reference_genome data/human_chr22/GRCh38.chr22.fa --sqanti_prefix data/human_chr22/rat_human_chr22 --n_random_isoforms 50 --output reference_data_chr22/human.chr22 ```

Step 2: generate expression profile based on real PacBio CCS data

```quantify.py --fastq data/human_chr22/Human.PacBio.ENCFF.chr22.fq -t 16 --reference_transcripts reference_data_chr22/human.chr22.transcripts.fasta --mandatory reference_data_chr22/human.chr22.novel_isoforms.tsv --output reference_data_chr22/human.chr22.counts.tsv```

Step 3: simulate data

```simulate.py --reference_prefix reference_data_chr22/human.chr22 --counts reference_data_chr22/human.chr22.counts.tsv -t 16 --output chr22_simulated```

The output files will be stored in `chr22_simulated` folder. Output description is provided in the following section.

## Output

- `simulator.log` log file
- `Illumina.simulated_1.fq`, `Illumina.simulated_2.fq` Illumina paired end reads
- `Illumina.simulated.sim.genes.results`, `Illumina.simulated.sim.isoforms.results` de facto expression values for Illumina reads
- `PacBio.simulated.fasta` simulated PacBio CCS reads
- `PacBio.simulated.isoform_counts.tsv` de facto counts for every isoform
- `PacBio.simulated.read_to_isoform.tsv` read ID to isofrom ID table
- `PacBio.simulated.tsv` internal IsoSeqSim file (detailed information on simulated isoforms)
-


## Reference data

Human ([gencode v37](https://www.gencodegenes.org/human/)):
- ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz
- ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.transcripts.fa.gz

Mouse ([gencode M26](https://www.gencodegenes.org/mouse/)):
- ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/gencode.vM26.annotation.gtf.gz
- ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/gencode.vM26.transcripts.fa.gz
