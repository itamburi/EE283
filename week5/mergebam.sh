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

# generated a prefix file with this command:
# ls *.bam | sed -E 's/_([0-9]+)\.sort\.bam//' | sort -u > prefix.txt

prefix=`head -n $SLURM_ARRAY_TASK_ID ${bamdir}/prefix.txt | tail -n 1`

# use printf to enumerate the replicates, and pass them into the argument of files to merge in samtools merge
reps=$(printf '%s ' ${bamdir}/${prefix}*.sort.bam)

samtools merge -u - ${reps} | samtools sort -o ${out}/${prefix}.mrg.bam


java -jar /opt/apps/picard-tools/2.27.1/picard.jar AddOrReplaceReadGroups\
	I=${out}/${prefix}.mrg.bam \
	O=${out}/${prefix}.RG.bam \
	SORT_ORDER=coordinate \
	RGPL=illumina \
	RGPU=D109LACXX \
	RGLB=Lib1 \
	RGID=${prefix} \
	RGSM=${prefix} \
	VALIDATION_STRINGENCY=LENIENT

java -jar /opt/apps/picard-tools/2.27.1/picard.jar MarkDuplicates REMOVE_DUPLICATES=true\
	I=${out}/${prefix}.RG.bam \
	O=${out}/${prefix}.dedup.bam \
	M=${out}/${prefix}_mrkdup_metrics.txt

samtools index ${out}/${prefix}.dedup.bam

ref="/data/class/ecoevo283/public/Bioinformatics_Course/ref/dmel-all-chromosome-r6.13.fasta"

java -jar /opt/apps/gatk/4.2.6.1/gatk HaplotypeCaller -R $ref \
	-I ${out}/${prefix}.dedup.bam \
	--minimum-mapping-quality 30 \
	-ERC GVCF \
	-O ${vcf}/${prefix}.g.vcf.gz

java -jar /opt/apps/gatk/4.2.6.1/gatk CombineGVCFs -R $ref $(printf -- '-V %s ' ${vcf}/*.g.vcf.gz) -O ${vcf}/allsample.g.vcf.gz





