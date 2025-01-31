#!/bin/bash
#SBATCH --job-name=bigwig    ## Name of the job.
#SBATCH -A class-ecoevo283 ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --cpus-per-task=16  ## number of cores the job needs
#SBATCH --mem=16GB

# source my conda software location

module load miniconda3/24.9.2
source /opt/apps/miniconda3/24.9.2/etc/profile.d/conda.sh
conda activate deeptools

# dirs for indexed bam files
A4="/dfs6/pub/itamburi/ee283/week3/prob4/DNAbam/ADL06_1.sort.bam"
A5="/dfs6/pub/itamburi/ee283/week3/prob4/DNAbam/ADL09_1.sort.bam"
dir="/dfs6/pub/itamburi/ee283/week4"

# prob3
bamCoverage -b $A4 -o ${dir}/A4/A4.bedgraph --normalizeUsing RPKM --region X:1903000:1905000 --binSize 10 --outFileFormat bedgraph
bamCoverage -b $A5 -o ${dir}/A5/A5.bedgraph --normalizeUsing RPKM --region X:1903000:1905000 --binSize 10 --outFileFormat bedgraph

# prob4
bamCoverage -b $A4 -o ${dir}/A4/A4_ext.bedgraph --extendReads 150 --normalizeUsing RPKM --region X:1903000:1905000 --binSize 10 --outFileFormat bedgraph
bamCoverage -b $A5 -o ${dir}/A5/A5_ext.bedgraph --extendReads 150 --normalizeUsing RPKM --region X:1903000:1905000 --binSize 10 --outFileFormat bedgraph
