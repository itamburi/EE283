#!/bin/sh
#SBATCH -A CLASS-ECOEVO283
#SBATCH --job-name=featurecounts    # Job name 
#SBATCH --cpus-per-task 64
#SBATCH --mem-per-cpu=3G
#SBATCH --output=/dfs6/pub/itamburi/ee283/week7/log/%x_%j.out   # log based on $
#SBATCH --error=/dfs6/pub/itamburi/ee283/week7/log/%x_%j.err    # error based o$

module load subread/2.0.3

gtf="/dfs6/pub/itamburi/ee283/datalinks/refs/dm6.gtf"
dir="/dfs6/pub/itamburi/ee283/week7"
# the gtf should match the genome you aligned to
# coordinates and chromosome names

# the program expects a space delimited set of files...
myfile=`cat ${dir}/shortRNAseq.names.txt | tr "\n" " "`
featureCounts -p -T 8 -t exon -g gene_id -Q 30 -F GTF -a $gtf -o ${dir}/fly_counts.txt $myfile
