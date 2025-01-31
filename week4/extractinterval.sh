module load samtools/1.15.1
module load bedtools2/2.30.0

# dirs for indexed bam files
A4="/dfs6/pub/itamburi/ee283/week3/prob4/DNAbam/ADL06_1.sort.bam"
A5="/dfs6/pub/itamburi/ee283/week3/prob4/DNAbam/ADL09_1.sort.bam"
dir="/dfs6/pub/itamburi/ee283/week4"

# *** Prob 1: extract all reads from an interval in a genome from A4 and A5 Bam files ****

# samtools idxstats [bamfile] to determine chromosome formatting 
interval="X:1880000-2000000"

# extract the read IDs that we will search for

samtools view $A4 $interval | cut -f1 > ${dir}/A4/A4.IDs.txt
samtools view $A5 $interval | cut -f1 > ${dir}/A5/A5.IDs.txt


# now we want to pull out these specific read IDs and save them in a .fa file
# since a BAM file contains one read, and its ID and metrics, per line, we can grep for the IDs in the text file
# then we use awk to extract only specific cols
# name=col1, sequence=col10, chromosome=col3
# note, we are excluding matches to chromosome "X" since this is where our interval is from

# All reads that map to chromosome X
samtools view $A4 | grep -f ${dir}/A4/A4.IDs.txt |\
    awk '{if($3 == "X"){printf(">%s\n%s\n",$1,$10)}}' >${dir}/A4/A4_X.fa

samtools view $A5 | grep -f ${dir}/A5/A5.IDs.txt |\
    awk '{if($3 == "X"){printf(">%s\n%s\n",$1,$10)}}' >${dir}/A5/A5_X.fa


# All reads that map to elsewhere besides chromX
samtools view $A4 | grep -f ${dir}/A4/A4.IDs.txt |\
    awk '{if($3 != "X"){printf(">%s\n%s\n",$1,$10)}}' >${dir}/A4/A4_other.fa

samtools view $A5 | grep -f ${dir}/A5/A5.IDs.txt |\
    awk '{if($3 != "X"){printf(">%s\n%s\n",$1,$10)}}' >${dir}/A5/A5_other.fa


#cat A4.fa | grep ">" | wc -l
#cat A5.fa | grep ">" | wc -l
#wc -l A4.IDs.txt
#wc -l A5.IDs.txt

# conda activate  spades before running
# need to specify full path to spades install
python /data/homezvol0/itamburi/.conda/envs/spades/bin/spades.py -o "/dfs6/pub/itamburi/ee283/week4/A4/assembly" -s ${dir}/A4/A4_other.fa --isolate > ${dir}/A4/assembly/A4.messages.txt
python /data/homezvol0/itamburi/.conda/envs/spades/bin/spades.py -o "/dfs6/pub/itamburi/ee283/week4/A5/assembly" -s ${dir}/A5/A5_other.fa --isolate > ${dir}/A5/assembly/A5.messages.txt

#cat spades_sssembly/contigs.fasta
