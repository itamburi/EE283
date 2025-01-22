module load fastqc/0.11.9 
module load trimmomatic/0.39
module load java/1.8.0

# *** Take a pair (i.e., READ1 and READ2) gDNA paired end sequence files and run fastqc on it
dir="/dfs6/pub/itamburi/ee283"

READ1=${dir}/datalinks/DNAseq/ADL06_1_1.fq.gz
READ2=${dir}/datalinks/DNAseq/ADL06_1_2.fq.gz

fastqc ${READ1} ${READ2} -o ${dir}/week3/fastqc1_out


# *** Run trimmomatic and rerun fastqc

# Trimmomatic
dirout="/dfs6/pub/itamburi/ee283/week3/prob1/trimmomatic_out"
trimmodir="/opt/apps/trimmomatic/0.39"
java -jar ${trimmodir}/trimmomatic-0.39.jar PE -threads 4 \
    ${READ1} ${READ2} \
    ${dirout}/ALD06_R1_paired.fq ${dirout}/ALD06_R1_unpaired.fq \
    ${dirout}/ALD06_R2_paired.fq ${dirout}/ALD06_R2_unpaired.fq \
    ILLUMINACLIP:${trimmodir}/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Rerun fastqc
fastqc ${dirout}/ALD06_R1_paired.fq ${dirout}/ALD06_R2_paired.fq -o ${dir}/week3/prob1/fastqc2_out


# *** Describe anything interesting

# Examining each .html report file produced by the first fastqc run, the per base sequence quality seems to decline at around a read position of 30-40 bp 

# Examining the .html report file produced by the fastqc on the reads processed by trimmomatic, the per base sequence quality remains high (>30) across the position in read (bp)

# Note: the trimmed reads files produced by trimmomatic would be contained in the trimmomatic_out folder, but I have added these to .gitignore because they exceed the github upload limit of 100mb
