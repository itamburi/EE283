#!/bin/sh
#SBATCH -A CLASS-ECOEVO283
#SBATCH --job-name=index.dm6    # Job name                  
#SBATCH --cpus-per-task 8
#SBATCH --mem=16gb
#SBATCH --output=/dfs6/pub/itamburi/ee283/week3/prob3/prob3-%J.log   # Output and error log 


module load picard-tools/2.27.1
module load bwa/0.7.8
module load samtools/1.15.1  

# load other modules
ref="/dfs6/pub/itamburi/ee283/datalinks/refs/dm6.fasta"

bwa index $ref -p "/dfs6/pub/itamburi/ee283/datalinks/refs/bwaidx"
samtools faidx $ref -o "/dfs6/pub/itamburi/ee283/datalinks/refs/samidx.gz"

# what is wrong with the command below?
java -jar /opt/apps/picard-tools/2.27.1/picard.jar \
	CreateSequenceDictionary R=$ref O="/dfs6/pub/itamburi/ee283/datalinks/refs/dm6.dict"









