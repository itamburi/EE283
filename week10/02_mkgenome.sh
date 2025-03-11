#!/bin/sh
#SBATCH -A mseldin_lab
#SBATCH --job-name=mkgenome    # Job name                  
#SBATCH --cpus-per-task 32
#SBATCH --mem=32G                    
#SBATCH --output=%J.log   # Output and error log 


ref='/dfs6/pub/itamburi/ee283/week10/susscrofa_ncbi/GCF_000003025.6'
gdir='/dfs6/pub/itamburi/ee283/week10/genomeidx'

module load star/2.7.10a
STAR --runMode genomeGenerate \
	--genomeDir ${gdir} \
	--genomeFastaFiles ${ref}/GCF_000003025.6_Sscrofa11.1_genomic.fna \
	--sjdbGTFfile ${ref}/genomic.gff \
	--runThreadN 5 &> log.txt
