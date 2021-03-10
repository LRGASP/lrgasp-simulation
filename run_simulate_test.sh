mkdir -p similation_test

tar -xzf data/human_chr22.tar.gz -C similation_test/

python prepare_reference_data.py \
  --reference_annotation similation_test/human_chr22/gencode.v36.annotation.chr22.gtf \
  --reference_transcripts similation_test/human_chr22/gencode.v36.transcripts.chr22.fa \
  --reference_genome similation_test/human_chr22/GRCh38.chr22.fa \
  --sqanti_prefix similation_test/human_chr22/rat_human_chr22 \
  --n_random_isoforms 50 \
  --output similation_test/reference_data_chr22/human.chr22

python quantify.py \
  --fastq similation_test/human_chr22/Human.PacBio.ENCFF.chr22.fq \
  -t 16 \
  --reference_transcripts similation_test/reference_data_chr22/human.chr22.transcripts.fasta \
  --mandatory similation_test/reference_data_chr22/human.chr22.novel_isoforms.tsv \
  --output similation_test/reference_data_chr22/human.chr22.counts.tsv

python simulate.py \
  --reference_prefix similation_test/reference_data_chr22/human.chr22 \
  --counts similation_test/reference_data_chr22/human.chr22.counts.tsv \
  -t 16 --test_mode \
  --output similation_test/chr22_simulated/