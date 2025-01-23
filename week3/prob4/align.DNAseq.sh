#!/bin/sh
#SBATCH -A CLASS-ECOEVO283
#SBATCH --job-name=aligndna    # Job name                  
#SBATCH --cpus-per-task 8
#SBATCH --array=1-12  
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/dfs6/pub/itamburi/ee283/week3/prob4/DNAlog/%J.log   # Output and error log 

module load bwa/0.7.8
module load samtools/1.15.1  

# bwaidx directtory contains index files with prefix bwaidx.* ...  we need to address up to the prefix for bwa
idx="/dfs6/pub/itamburi/ee283/datalinks/refs/bwaidx/bwaidx"
data="/dfs6/pub/itamburi/ee283/datalinks/DNAseq"
out="/dfs6/pub/itamburi/ee283/week3/prob4/DNAbam"

# generated prefixes.txt by running the following within DNAseq directory
# ls *.fq.gz | sed -E 's/_[^_]+\.fq\.gz$//' | sort -u > prefixes.txt

prefix=`head -n $SLURM_ARRAY_TASK_ID ${data}/prefixes.txt | tail -n 1`


bwa mem -t $SLURM_CPUS_PER_TASK -M $idx ${data}/${prefix}_1.fq.gz ${data}/${prefix}_2.fq.gz \
	| samtools sort - -@ $SLURM_CPUS_PER_TASK -m 4G -o ${out}/${prefix}.sort.bam












