#!/bin/sh
#SBATCH -A CLASS-ECOEVO283
#SBATCH --cpus-per-task 12
#SBATCH --mem-per-cpu=3G
#SBATCH --job-name=biallelic
#SBATCH --output=/dfs6/pub/itamburi/ee283/week5/log/%x.out   # log based on jobname
#SBATCH --error=/dfs6/pub/itamburi/ee283/week5/log/%x.err    # error based on jobname

module load vcftools/0.1.16

vcf="/dfs6/pub/itamburi/ee283/week5/vcf"
dir="/dfs6/pub/itamburi/ee283/week5"

vcftools --gzvcf ${vcf}/filt3_variants.vcf.gz \
	--chr X \
	--from-bp 1 \
	--to-bp 1000000 \
	--max-alleles 2 \
	--min-alleles 2 \
	--remove-indels \
	--012 \
	--out ${dir}/X_1Mb_snps
