#!/bin/sh
#SBATCH -A CLASS-ECOEVO283
#SBATCH --job-name=aligndna    # Job name                  
#SBATCH --cpus-per-task 8
#SBATCH --array=1-12  
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/dfs6/pub/itamburi/ee283/week3/prob4/DNAlog/samtoolsindex_log/%J.log   # Output and error log 

module load samtools/1.15.1  

out="/dfs6/pub/itamburi/ee283/week3/prob4/DNAbam"
data="/dfs6/pub/itamburi/ee283/datalinks/DNAseq"
prefix=`head -n $SLURM_ARRAY_TASK_ID ${data}/prefixes.txt | tail -n 1`

samtools index ${out}/${prefix}.sort.bam 











