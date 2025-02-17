#!/bin/sh
#SBATCH -A CLASS-ECOEVO283
#SBATCH --job-name=call_peaks    # Job name 
#SBATCH --cpus-per-task 8
#SBATCH --array=1-4 
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/dfs6/pub/itamburi/ee283/week6/log/%x_%j.out   # log based on jobname
#SBATCH --error=/dfs6/pub/itamburi/ee283/week6/log/%x_%j.err    # error based on jobname


module load samtools/1.15.1
module load bedtools2/2.30.0
module load ucsc-tools/v429

wd="/dfs6/pub/itamburi/ee283/week6"
samples="/dfs6/pub/itamburi/ee283/week6/atac.samples.txt"

prefix=`head -n $SLURM_ARRAY_TASK_ID ${samples} | tail -n 1`


samtools merge -o condition1.bam sample1.bam sample2.bam sample3.bam

bedtools bamtobed -i condition1.bam | \
	awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' \
	> condition1.tn5.bed

macs2 callpeak -t condition1.tn5.bed \
	-n condition1 \
	-f BED -g mm -q 0.01 --nomodel \
	--shift -75 --extsize 150 \
	--call-summits --keep-dup all -B --broad

# I cd into my scripts folder and download a function ucsctools doesn't have.
# Run in interactive session:
# rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64.v369/fetchChromSizes 

${wd}/scripts/fetchChromSizes ${ref}/dm6 > ${ref}/dm6.chrom.sizes

LC_COLLATE=C sort -k1,1 -k2,2n ${wd}/bdg.files/broad_treat_pileup.bdg > ${wd}/bdg.files/broad_treat_pileup.sorted.bdg

awk 'NR==FNR {chr_size[$1]=$2; next} $3 > chr_size[$1]' ${ref}/dm6.chrom.sizes ${wd}/bdg.files/broad_treat_pileup.sorted.bdg > ${wd}/bdg.files/broad_treat_pileup.safeends.bdg

bedGraphToBigWig broad_treat_pileup.safeends.bdg ref/dm6.chrom_sizes broad_peaks.bw



###########
module load java/1.8.0
module load samtools/1.15.1
module load picard-tools/2.27.1

samples="/dfs6/pub/itamburi/ee283/week6/atac.samples.txt"
atacbam="/dfs6/pub/itamburi/ee283/week3/prob4/ATACbam"
bamout="/dfs6/pub/itamburi/ee283/week6/bam.processed"

prefix=`head -n $SLURM_ARRAY_TASK_ID ${samples} | tail -n 1`

# need indexed file to view chr X in the next step
samtools index -c ${atacbam}/${prefix}.sort.bam

samtools view -q 30 -b ${atacbam}/${prefix}.sort.bam X | \
	samtools sort -O BAM -o ${bamout}/${prefix}.chrX.sort.bam #probably dont need to sort again, but just in case...

java -jar /opt/apps/picard-tools/2.27.1/picard.jar MarkDuplicates \
	I=${bamout}/${prefix}.chrX.sort.bam \
	O=${bamout}/${prefix}.chrX.dd.bam \
	M=${bamout}/${prefix}.chrX.dd.metrics \
	REMOVE_DUPLICATES=true




