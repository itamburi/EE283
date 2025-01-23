#!/bin/sh
#SBATCH -A CLASS-ECOEVO283
#SBATCH --job-name=alignrna    # Job name                  
#SBATCH --cpus-per-task 8
#SBATCH --mem=16G
#SBATCH --output=/dfs6/pub/itamburi/ee283/week3/prob5/indexlog-%J.log   # Output and error log 

module load hisat2/2.2.1
dir="/dfs6/pub/itamburi/ee283/datalinks/refs"
fa="${dir}/dm6.fasta"
gtf="${dir}/dm6.gtf"

python hisat2_extract_splice_sites.py $gtf > "${dir}/hisat/dm6.ss"
python hisat2_extract_exons.py $gtf > "${dir}/hisat/dm6.exon"
hisat2-build -p 8 --exon ${dir}/hisat/dm6.exon --ss ${dir}/hisat/dm6.ss $fa ${dir}/hisat/dm6_trans











