#!/bin/sh
#SBATCH -A CLASS-ECOEVO283
#SBATCH --job-name=call_peaks    # Job name 
#SBATCH --cpus-per-task 32
#SBATCH --mem-per-cpu=3G
#SBATCH --output=/dfs6/pub/itamburi/ee283/week6/log/%x_%j.out   # log based on jobname
#SBATCH --error=/dfs6/pub/itamburi/ee283/week6/log/%x_%j.err    # error based on jobname


module load samtools/1.15.1
module load bedtools2/2.30.0
module load ucsc-tools/v429

module load miniconda3/24.9.2
source /opt/apps/miniconda3/24.9.2/etc/profile.d/conda.sh
conda activate macs2

ref="/dfs6/pub/itamburi/ee283/datalinks/refs"
wd="/dfs6/pub/itamburi/ee283/week6"
bam="/dfs6/pub/itamburi/ee283/week6/bam.processed"
out="/dfs6/pub/itamburi/ee283/week6/pkfiles"
samples="/dfs6/pub/itamburi/ee283/week6/atac.samples.txt"


samtools merge -f -o ${out}/A4ED_merged.bam \
	${bam}/P004_A4_ED_2.chrX.dd.bam \
	${bam}/P005_A4_ED_3.chrX.dd.bam \
	${bam}/P006_A4_ED_4.chrX.dd.bam

bedtools bamtobed -i ${out}/A4ED_merged.bam | \
	awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' \
	> ${out}/A4ED_merged.tn5.bed

macs2 callpeak -t ${out}/A4ED_merged.tn5.bed \
	--outdir ${out} \
	-n A4ED \
	-f BED -g mm -q 0.01 --nomodel \
	--shift -75 --extsize 150 \
	--broad --keep-dup all --bdg



# I cd into my scripts folder and download a function ucsctools doesn't have.
# Run in interactive session:
# rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64.v369/fetchChromSizes . 
# ${wd}/scripts/fetchChromSizes ${ref}/dm6.fasta > ${ref}/dm6.chrom.sizes
# ... this command didnt seem to work. So isntead I downloaded the chrom sizes file from UCSC and reformmated by removing the chr:
#wget -O - http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes | sed 's/^chr//' > dm6.chrom.sizes.ucsc

LC_COLLATE=C sort -k1,1 -k2,2n ${out}/A4ED_treat_pileup.bdg > ${out}/A4ED_treat_pileup.sorted.bdg

# awk 'NR==FNR {chr_size[$1]=$2; next} $3 > chr_size[$1]' ${ref}/dm6.chrom.sizes.ucsc ${out}/A4ED_treat_pileup.sorted.bdg > ${out}/A4ED_treat_pileup.safeends.bdg
# ... this didnt work, probably because the .sorted.bdg only contains X chromosome. instead use:
awk 'NR==FNR {if ($1 == "X") chr_size[$1]=$2; next} $1 == "X" && $3 <= chr_size[$1]' ${ref}/dm6.chrom.sizes.ucsc ${out}/A4ED_treat_pileup.sorted.bdg > ${out}/A4ED_treat_pileup.safeends.bdg 

#awk '{print "chr" $1 "\t" $2}' ${ref}/dm6.chrom.sizes.ucsc > ${ref}/dm6_with_chr.chrom.sizes.ucsc
#awk '{if ($1 ~ /^[0-9]+$/) $1 = "chr" $1; else if ($1 !~ /^chr/) $1 = "chr" $1; print $0}' ${out}/A4ED_treat_pileup.safeends.bdg > ${out}/A4ED_treat_pileup_with_chr.bdg
# this one tab delimits: awk 'BEGIN {OFS="\t"} {if ($1 ~ /^[0-9]+$/) $1 = "chr" $1; else if ($1 !~ /^chr/) $1 = "chr" $1; print $0}' ${out}/A4ED_treat_pileup.safeends.bdg > ${out}/A4ED_treat_pileup_with_chr.bdg
# .. but it doesnt seem to have any change on the resulting .bw below

bedGraphToBigWig ${out}/A4ED_treat_pileup_with_chr.bdg ${ref}/dm6_with_chr.chrom.sizes.ucsc ${out}/A4ED_broad_peaks_with_chr.bw
#bedGraphToBigWig ${out}/A4ED_treat_pileup.safeends.bdg ${ref}/dm6.chrom.sizes.ucsc ${out}/A4ED_broad_peaks.bw



