#!/bin/sh
#SBATCH -A CLASS-ECOEVO283
#SBATCH --job-name=qc_atac    # Job name 
#SBATCH --cpus-per-task 32
#SBATCH --array=1-6 
#SBATCH --mem-per-cpu=3G
#SBATCH --output=/dfs6/pub/itamburi/ee283/week6/log/%x_%j.out   # log based on jobname
#SBATCH --error=/dfs6/pub/itamburi/ee283/week6/log/%x_%j.err    # error based on jobname

# We will run QC on only the X Chromosome. And we will only check out a couple of samples by arraying 1-6  for expediency

module load java/1.8.0
module load samtools/1.15.1
module load picard-tools/2.27.1

samples="/dfs6/pub/itamburi/ee283/week6/atac.samples.txt"
atacbam="/dfs6/pub/itamburi/ee283/week3/prob4/ATACbam"
bamout="/dfs6/pub/itamburi/ee283/week6/bam.processed"

prefix=`head -n $SLURM_ARRAY_TASK_ID ${samples} | tail -n 1`

# need indexed file to view chr X in the next step
samtools index -c ${atacbam}/${prefix}.sort.bam

samtools view -q 30 -b ${atacbam}/${prefix}.sort.bam X | \
	samtools sort -O BAM -o ${bamout}/${prefix}.chrX.sort.bam #probably dont need to sort again, but just in case...

java -jar /opt/apps/picard-tools/2.27.1/picard.jar MarkDuplicates \
	I=${bamout}/${prefix}.chrX.sort.bam \
	O=${bamout}/${prefix}.chrX.dd.bam \
	M=${bamout}/${prefix}.chrX.dd.metrics \
	REMOVE_DUPLICATES=true


samtools index ${bamout}/${prefix}.chrX.dd.bam

rm ${bamout}/${prefix}.chrX.sort.bam

