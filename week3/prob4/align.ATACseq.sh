#!/bin/sh
#SBATCH -A CLASS-ECOEVO283
#SBATCH --job-name=alignATAC    # Job name                  
#SBATCH --cpus-per-task 4
#SBATCH --array=1-24  
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/dfs6/pub/itamburi/ee283/week3/prob4/ATAClog/%J.log   # Output and error log 

module load bwa/0.7.8
module load samtools/1.15.1  

idx="/dfs6/pub/itamburi/ee283/datalinks/refs/bwaidx/bwaidx"
data="/dfs6/pub/itamburi/ee283/datalinks/ATACseq"
bamout="/dfs6/pub/itamburi/ee283/week3/prob4/ATACbam"

# generated prefixes.txt by running the following within DNAseq directory
# ls *.fq.gz | sed -E 's/_[^_]+\.fq\.gz$//' | sort -u > prefixes.txt
prefix=`head -n $SLURM_ARRAY_TASK_ID ${data}/prefixes.txt | tail -n 1`


# *** Run trimmomatic
trimout="/dfs6/pub/itamburi/ee283/week3/prob4/ATAC_trimmo"
loc="/opt/apps/trimmomatic/0.39"
java -jar ${loc}/trimmomatic-0.39.jar PE -threads 4 \
    ${data}/${prefix}_R1.fq.gz ${data}/${prefix}_R2.fq.gz \
    ${trimout}/${prefix}_R1_paired.fq ${trimout}/${prefix}_R1_unpaired.fq \
    ${trimout}/${prefix}_R2_paired.fq ${trimout}/${prefix}_R2_unpaired.fq \
    ILLUMINACLIP:${loc}/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# *** Align with BWA
bwa mem -t $SLURM_CPUS_PER_TASK -M $idx ${trimout}/${prefix}_R1_paired.fq ${trimout}/${prefix}_R2_paired.fq \
	| samtools sort - -@ $SLURM_CPUS_PER_TASK -m 4G -o ${bamout}/${prefix}.sort.bam












