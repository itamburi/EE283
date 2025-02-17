#!/bin/sh
#SBATCH -A CLASS-ECOEVO283
#SBATCH --job-name=mergevcf    # Job name                  
#SBATCH --cpus-per-task 12
#SBATCH --mem-per-cpu=3G
#SBATCH --output=/dfs6/pub/itamburi/ee283/week5/log/%J.log   # Output and error log 

module load  bcftools/1.15.1

vcf="/dfs6/pub/itamburi/ee283/week5/vcf"

bcftools filter -i 'FS<40.0 && SOR<3 && MQ>40.0 && MQRankSum>-5.0 && MQRankSum<5 && QD>2.0 && ReadPosRankSum>-4.0 && INFO/DP<16000' -O z -o ${vcf}/filt1_variants.vcf.gz ${vcf}/all_variants.vcf.gz

bcftools filter -S . -e 'FMT/DP<3 | FMT/GQ<20' -O z -o ${vcf}/filt2_variants.vcf.gz ${vcf}/filt1_variants.vcf.gz

bcftools filter -e 'AC==0 || AC==AN' --SnpGap 10 ${vcf}/filt2_variants.vcf.gz | \
	bcftools view -m2 -M2 -v snps -O z -o ${vcf}/filt3_variants.vcf.gz
