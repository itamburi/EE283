#!/bin/sh
#SBATCH -A CLASS-ECOEVO283
#SBATCH --job-name=combinevcf    # Job name                  
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/dfs6/pub/itamburi/ee283/week5/log/%J.log   # Output and error log 

module load java/1.8.0
module load gatk/4.2.6.1
module load picard-tools/2.27.1
module load samtools/1.15.1

ref="/dfs6/pub/itamburi/ee283/datalinks/refs/dm6.fasta"
vcf="/dfs6/pub/itamburi/ee283/week5/vcf"
/opt/apps/gatk/4.2.6.1/gatk CombineGVCFs -R $ref $(printf -- '-V %s ' ${vcf}/*.g.vcf.gz) -O ${vcf}/allsample.g.vcf.gz
# the printif line expands all the files in the vcf dir to -V [filename] as inputs to combine







