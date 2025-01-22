#### problem2.sh ####

module load fastqc/0.11.9 
module load trimmomatic/0.39
module load java/1.8.0

# Run fastqc on a pair of ATACseq fastq files
dir="/dfs6/pub/itamburi/ee283"

READ1=${dir}/datalinks/ATACseq/P004_A4_ED_2_R1.fq.gz
READ2=${dir}/datalinks/ATACseq/P004_A4_ED_2_R2.fq.gz

# NOTE: initially my symlinks didnt have .fq.gz extensions, and fastqc couldnt parse them correctly. Updated my symlinks to have the .fq.gz extension
fastqc ${READ1} ${READ2} -o ${dir}/week3/prob2/fastqc1_out


# *** Run trimmomatic
# Trimmomatic
dirout="/dfs6/pub/itamburi/ee283/week3/prob2/trimmomatic_out"
trimmodir="/opt/apps/trimmomatic/0.39"
java -jar ${trimmodir}/trimmomatic-0.39.jar PE -threads 4 \
    ${READ1} ${READ2} \
    ${dirout}/P004_A4_ED_2_R1_paired.fq ${dirout}/P004_A4_ED_2_R1_unpaired.fq \
    ${dirout}/P004_A4_ED_2_R2_paired.fq ${dirout}/P004_A4_ED_2_R2_unpaired.fq \
    ILLUMINACLIP:${trimmodir}/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


# *** Rerun fastqc
fastqc ${dirout}/P004_A4_ED_2_R1_paired.fq ${dirout}/P004_A4_ED_2_R2_paired.fq -o ${dir}/week3/prob2/fastqc2_out


# *** Describe anything interesting

# once again, the qualities are imporved according to the .html reports
