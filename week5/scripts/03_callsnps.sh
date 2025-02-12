#!/bin/sh
#SBATCH -A CLASS-ECOEVO283
#SBATCH --job-name=GATK    # Job name                  
#SBATCH --cpus-per-task 32
#SBATCH --array=1-7
#SBATCH --mem-per-cpu=3G
#SBATCH --output=/dfs6/pub/itamburi/ee283/week5/log/%J.log   # Output and error log 

module load java/1.8.0
module load gatk/4.2.6.1 
module load bcftools/1.15.1
module load samtools/1.15.1 

ref="/dfs6/pub/itamburi/ee283/datalinks/refs/dm6.fasta"
dir="/dfs6/pub/itamburi/ee283/week5"
vcf="/dfs6/pub/itamburi/ee283/week5/vcf"

chr=`head -n $SLURM_ARRAY_TASK_ID ${dir}/chrome.names.txt | tail -n 1`

/opt/apps/gatk/4.1.9.0/gatk --java-options "-Xmx3g" GenotypeGVCFs \
	-R $ref \
	-V ${vcf}/allsample.g.vcf.gz \
	--intervals ${chr} \
	-stand-call-conf 5 \
	-O ${vcf}/result.${chr}.vcf



# Compress the VCF with BGZip
bgzip ${vcf}/result.${chr}.vcf

# Index the BGZip-compressed VCF
bcftools index ${vcf}/result.${chr}.vcf.gz

# resulting script produces .idx and .csi files. Removing the .idx, .csi is produced from bcftools index I  beleive when inputs are large

