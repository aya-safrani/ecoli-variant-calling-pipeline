#!/usr/bin/env bash
# Variant calling pipeline on simulated E. coli K-12 reads
# Run from the project root. Tool versions: see results/versions_outils.txt
# Note: the samtools steps (section 4) are rebuilt from the intermediate files kept in align/.
set -euo pipefail

THREADS=4
REF=data/ecoli_ref.fasta
# SnpEff database for E. coli K-12 MG1655 (check available names with: snpEff databases | grep -i mg1655)
DB="Escherichia_coli_str_k_12_substr_mg1655_gca_000005845"

mkdir -p qc/raw qc/trimmed trimmed align variants annotation results

# 0. Environment
# conda create -n tp_genomique --solver=libmamba --override-channels -c conda-forge -c bioconda \
#   fastqc multiqc trimmomatic bwa samtools bcftools snpeff wgsim -y

# 1. Simulate paired-end reads (460,000 pairs, 2x150 bp) with known variants
( cd data && wgsim -N 460000 -1 150 -2 150 -r 0.001 -e 0.001 \
    ecoli_ref.fasta sample_R1.fastq sample_R2.fastq > wgsim_truth_variants.txt )

# 2. Quality control of raw reads
fastqc -t "$THREADS" data/sample_R1.fastq data/sample_R2.fastq -o qc/raw/
multiqc qc/raw/ -o qc/raw/

# 3. Trimming
trimmomatic PE -threads "$THREADS" -phred33 \
  data/sample_R1.fastq data/sample_R2.fastq \
  trimmed/R1_paired.fastq trimmed/R1_unpaired.fastq \
  trimmed/R2_paired.fastq trimmed/R2_unpaired.fastq \
  LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50 \
  2>&1 | tee results/q2_trimmomatic.log

# 4. Mapping and BAM processing
bwa index "$REF"
bwa mem -t "$THREADS" -R '@RG\tID:sample1\tSM:sample1\tPL:ILLUMINA' \
  "$REF" trimmed/R1_paired.fastq trimmed/R2_paired.fastq > align/sample.sam
samtools view -b align/sample.sam > align/sample.bam
samtools sort -n -o align/sample.namesort.bam align/sample.bam
samtools fixmate -m align/sample.namesort.bam align/sample.fixmate.bam
samtools sort -o align/sample.sorted.bam align/sample.fixmate.bam
samtools markdup align/sample.sorted.bam align/sample.dedup.bam
samtools index align/sample.dedup.bam
samtools flagstat align/sample.dedup.bam > results/q3_flagstat.txt

# 5. Variant calling (diploid model) and filtering
bcftools mpileup -f "$REF" align/sample.dedup.bam \
  | bcftools call -mv -Oz -o variants/sample.raw.vcf.gz
bcftools filter -e 'QUAL<30 || INFO/DP<10' \
  variants/sample.raw.vcf.gz -Oz -o variants/sample.filtered.vcf.gz
bcftools stats variants/sample.filtered.vcf.gz > variants/sample.stats.txt

# 5b. Haploid model (--ploidy 1), for comparison
bcftools mpileup -f "$REF" align/sample.dedup.bam \
  | bcftools call -mv --ploidy 1 -Oz -o variants/sample.raw.haploid.vcf.gz
bcftools filter -e 'QUAL<30 || INFO/DP<10' variants/sample.raw.haploid.vcf.gz \
  -Oz -o variants/sample.filtered.haploid.vcf.gz

# 6. Functional annotation with SnpEff
# chr_to_snpeff.txt maps the reference chromosome name to the one used by the SnpEff database
snpEff download -v "$DB"
bcftools annotate --rename-chrs annotation/chr_to_snpeff.txt \
  variants/sample.filtered.vcf.gz \
  | snpEff ann -v -stats annotation/snpEff_summary.html "$DB" \
  > annotation/sample.annotated.snpeff_names.vcf 2> results/q5_snpeff.log

# 7. Visual inspection of variants: open align/sample.dedup.bam in IGV (see igv_snapshot.png)
