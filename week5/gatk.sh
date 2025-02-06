#!/bin/sh
#SBATCH -A CLASS-ECOEVO283
#SBATCH --job-name=GATK    # Job name                  
#SBATCH --cpus-per-task 4
#SBATCH --array=1-4 		# there are 384 samples but we only need to do a few for the hw 
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/dfs6/pub/itamburi/ee283/week5/log/%J.log   # Output and error log 

module load java/1.8.0
module load gatk/4.2.6.1 
module load picard-tools/2.27.1  
module load samtools/1.15.1


bamdir="/dfs6/pub/itamburi/ee283/week3/prob4/DNAbam"
out="/dfs6/pub/itamburi/ee283/week5/bamout"
vcf="/dfs6/pub/itamburi/ee283/week5/vcf"


prefix=`head -n $SLURM_ARRAY_TASK_ID ${bamdir}/prefix.txt | tail -n 1`


java -jar /opt/apps/gatk/4.2.6.1/gatk HaplotypeCaller -R $ref \
	-I ${out}/${prefix}.dedup.bam \
	--minimum-mapping-quality 30 \
	-ERC GVCF \
	-O ${vcf}/${prefix}.g.vcf.gz

java -jar /opt/apps/gatk/4.2.6.1/gatk CombineGVCFs -R $ref $(printf -- '-V %s ' ${vcf}/*.g.vcf.gz) -O ${vcf}/allsample.g.vcf.gz





