#!/bin/sh
#SBATCH -A CLASS-ECOEVO283
#SBATCH --job-name=RNAseqAlign    # Job name                  
#SBATCH --cpus-per-task 32
#SBATCH --array=1-20 		# there are 384 samples but we only need to do a few for the hw 
#SBATCH --mem-per-cpu=3G
#SBATCH --error=/dfs6/pub/itamburi/ee283/week3/prob5/RNAlog/RNAseqAlign.%A_%a.err
#SBATCH --output=/dfs6/pub/itamburi/ee283/week3/prob5/RNAlog/RNAseqAlign.%A_%a.out

module load samtools/1.15.1  
module load hisat2/2.2.1

idx="/dfs6/pub/itamburi/ee283/datalinks/refs/bwaidx/bwaidx"
data="/dfs6/pub/itamburi/ee283/datalinks/RNAseq"
bamout="/dfs6/pub/itamburi/ee283/week3/prob5/RNAbam"

# generated prefixes.txt by running the following within RNAseq directory
# ls *.fq.gz | sed -E 's/_[^_]+\.fq\.gz$//' | sort -u > prefixes.txt
# ls | grep -E '21148|21286|22162|21297|21029|22052|22031|21293|22378|22390' | cut -d'_' -f1 | grep -E 'E|B' | sort | uniq > prefixes.subset.txt
prefix=`head -n $SLURM_ARRAY_TASK_ID ${data}/prefixes.subset.txt | tail -n 1`


idx="/dfs6/pub/itamburi/ee283/datalinks/refs/hisat/dm6_trans"

hisat2 -p 2 -x ${idx} -1 ${data}/${prefix}_R1.fq.gz -2 ${data}/${prefix}_R2.fq.gz \
	| samtools sort - -@ $SLURM_CPUS_PER_TASK -m 4G -o ${bamout}/${prefix}.sort.bam

samtools index ${bamout}/${prefix}.sort.bam 





