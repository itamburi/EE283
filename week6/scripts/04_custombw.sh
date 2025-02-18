#!/bin/sh
#SBATCH -A CLASS-ECOEVO283
#SBATCH --job-name=mkbw    # Job name 
#SBATCH --cpus-per-task 32
#SBATCH --mem-per-cpu=3G
#SBATCH --output=/dfs6/pub/itamburi/ee283/week6/log/%x_%j.out   # log based on jobname
#SBATCH --error=/dfs6/pub/itamburi/ee283/week6/log/%x_%j.err    # error based on jobname



module load samtools/1.15.1  
module load bedtools2/2.30.0  
module load ucsc-tools/v429

out="/dfs6/pub/itamburi/ee283/week6/pkfiles"
ref="/dfs6/pub/itamburi/ee283/datalinks/refs"

samtools sort ${out}/A4ED_merged.bam -o ${out}/A4ED_sorted.bam
samtools index ${out}/A4ED_sorted.bam

# count the total number of reads
# and create a scaling constant for each bam
chromsizes=${ref}/dm6.chrom.sizes.ucsc
Nreads=`samtools view -c -q 30 -F 4 ${out}/A4ED_sorted.bam`
Scale=`echo "1.0/($Nreads/1000000)" | bc -l`
samtools view -b ${out}/A4ED_sorted.bam | genomeCoverageBed -ibam - -g ${ref}/dm6.fasta -bg -scale $Scale > ${out}/A4ED.coverage
# these are the so called "kent-tools", useful for converting between various file types
bedGraphToBigWig ${out}/A4ED.coverage $chromsizes ${out}/A4ED.bw
