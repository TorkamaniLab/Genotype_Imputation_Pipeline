# VCF_build_check_tool
This tool identifies the genome build version (hg18, GRCh37, NCBI36, etc) that was used to make the input vcf file, and tells you the correct chain file to be used for mapping the input to GRCh37 (the same build used by 1000 Genomes or HRC reference panels).

#Requirements
R 3.3.2 or more recent versions

# Installation
Extract all *.bi files:

      gunzip *.gz

# Usage

      Rscript check_vcf_build.R input.vcf > output.txt

# Output example
      "Reading first 100K markers from chr1 in input file: /stsi/raqueld/vcf/ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out_bed.vcf..."
      "Reading chr1_GRCh37-hg19.hgTable.bi from /gpfs/home/raqueld/check_vcf_build"
      "Identifying matches between /gpfs/home/raqueld/check_vcf_build/chr1_GRCh37-hg19.hgTable.bi and /stsi/raqueld/vcf/ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out_bed.vcf"
      "Exact pos/ref/alt matches found: /stsi/raqueld/vcf/ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out_bed.vcf versus /gpfs/home/raqueld/check_vcf_build/chr1_GRCh37-hg19.hgTable.bi = 25"
      "Reading chr1_GRCh38-hg38.hgTable.bi from /gpfs/home/raqueld/check_vcf_build"
      "Identifying matches between /gpfs/home/raqueld/check_vcf_build/chr1_GRCh38-hg38.hgTable.bi and /stsi/raqueld/vcf/ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out_bed.vcf"
      "Exact pos/ref/alt matches found: /stsi/raqueld/vcf/ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out_bed.vcf versus /gpfs/home/raqueld/check_vcf_build/chr1_GRCh38-hg38.hgTable.bi = 31"
      "Reading chr1_NCBI34-hg16.hgTable.bi from /gpfs/home/raqueld/check_vcf_build"
      "Identifying matches between /gpfs/home/raqueld/check_vcf_build/chr1_NCBI34-hg16.hgTable.bi and /stsi/raqueld/vcf/ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out_bed.vcf"
      "Exact pos/ref/alt matches found: /stsi/raqueld/vcf/ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out_bed.vcf versus /gpfs/home/raqueld/check_vcf_build/chr1_NCBI34-hg16.hgTable.bi = 21"
      "Reading chr1_NCBI35-hg17.hgTable.bi from /gpfs/home/raqueld/check_vcf_build"
      "Identifying matches between /gpfs/home/raqueld/check_vcf_build/chr1_NCBI35-hg17.hgTable.bi and /stsi/raqueld/vcf/ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out_bed.vcf"
      "Exact pos/ref/alt matches found: /stsi/raqueld/vcf/ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out_bed.vcf versus /gpfs/home/raqueld/check_vcf_build/chr1_NCBI35-hg17.hgTable.bi = 21"
      "Reading chr1_NCBI36-hg18.hgTable.bi from /gpfs/home/raqueld/check_vcf_build"
      "Identifying matches between /gpfs/home/raqueld/check_vcf_build/chr1_NCBI36-hg18.hgTable.bi and /stsi/raqueld/vcf/ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out_bed.vcf"
      "Exact pos/ref/alt matches found: /stsi/raqueld/vcf/ARIC_PLINK_flagged_chromosomal_abnormalities_zeroed_out_bed.vcf versus /gpfs/home/raqueld/check_vcf_build/chr1_NCBI36-hg18.hgTable.bi = 15538"
      "Dataset is based on build: NCBI36/hg18. Exact pos/ref/alt matches: 15538. Total accidental matches from other builds, summed together: 98"
      "Use chain file(s): NCBI36_to_GRCh37.chain"
