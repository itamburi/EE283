#!/bin/sh
#SBATCH -A CLASS-ECOEVO283
#SBATCH --job-name=mergevcf    # Job name                  
#SBATCH --cpus-per-task 12
#SBATCH --mem-per-cpu=3G
#SBATCH --output=/dfs6/pub/itamburi/ee283/week5/log/%J.log   # Output and error log 

module load java/1.8.0
module load bcftools/1.15.1
module load picard-tools/2.27.1  

vcf="/dfs6/pub/itamburi/ee283/week5/vcf"

java -jar /opt/apps/picard-tools/2.27.1/picard.jar MergeVcfs $(printf 'I=%s ' ${vcf}/result.*.vcf.gz) O=${vcf}/all_variants.vcf.gz

bcftools index ${vcf}/all_variants.vcf.gz
