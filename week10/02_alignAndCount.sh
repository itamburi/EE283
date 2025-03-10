#!/bin/sh
#SBATCH -A mseldin_lab
#SBATCH --job-name=test    # Job name                  
#SBATCH --cpus-per-task 32
#SBATCH --array=1-10
#SBATCH --mem=3gb
#SBATCH --constraint=nvme
#SBATCH --output=/dfs6/pub/itamburi/star/test8/slurm/ldl-%J.log   # Output and error log 

# source my conda software location
source /data/homezvol0/itamburi/.bashrc
# activate conda enviornment with htseq
conda activate htseq

# load modules
module load star/2.7.10a
module load samtools/1.15.1
module load picard-tools/2.27.1

# set up enviornment variables
dir="/dfs6/pub/itamburi/ee283/week10"
prefix=`cat ${dir}/prefixes.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
reads=$(printf '%s ' ${dir}/datalinks/${prefix}_[12].fq.gz)

cd ${TMPDIR}
STAR --genomeDir ${dir}/genome \
        --readFilesIn ${reads} \
        --outSAMtype BAM Unsorted \
        --readFilesCommand zcat \
        --outFileNamePrefix ${prefix}_ \
        --outTmpDir ${TMPDIR}/${SLURM_ARRAY_TASK_ID} \
        --runThreadN 10 \

# Aligned.out.bam is the output from STAR based on outSAMtype
samtools sort ${prefix}_Aligned.out.bam -o ${TMPDIR}/${prefix}_sorted.bam

# markdupes and build counts matrix
picard MarkDuplicates \
    -INPUT ${TMPDIR}/${prefix}_sorted.bam \
    -OUTPUT ${TMPDIR}/${prefix}_dd.bam \
    -COMPRESSION_LEVEL 0 \
    -ASSUME_SORTED true \
    -SORTING_COLLECTION_SIZE_RATIO 0.8 \
    -METRICS_FILE ${TMPDIR}/${prefix}_mrkdp_met.txt \
    -VALIDATION_STRINGENCY SILENT \
    -REMOVE_DUPLICATES true


ref='/dfs6/pub/itamburi/ee283/week10/genome.ncbi/genomic.gtf'

htseq-count -f bam -r pos -i gene_id --stranded=no \
        ${TMPDIR}/${sample}_dd.bam \
        ${ref} \
> ${sample}_counts.txt


# move files back to CRSP
dir1="${dir}/counts"

if [ ! -d "$dir1" ]; then
    mkdir -p $dir1
    echo "Directory created: $dir1"
else
    echo "Directory already exists: $dir1"
fi

dir2="${dir}/starout"

if [ ! -d "$dir2" ]; then
    mkdir -p $dir2
    echo "Directory created: $dir2"
else
    echo "Directory already exists: $dir2"
fi

dir3="${dir}/sj"
if [ ! -d "$dir3" ]; then
    mkdir -p $dir3
    echo "Directory created: $dir3"
else
    echo "Directory already exists: $dir3"
fi

dir4="${dir}/mrkdupes"
if [ ! -d "$dir4" ]; then
    mkdir -p $dir4
    echo "Directory created: $dir4"
else
    echo "Directory already exists: $dir4"
fi


mv ${prefix}_counts.txt $dir1

mv ${prefix}_Log.out ${prefix}_Log.final.out ${prefix}_Log.progress.out $dir2
#mv ${sample}_Aligned.out.bam $dir2

mv ${prefix}_SJ.out.tab $dir3

mv ${prefix}_mrkdp_met.txt $dir4
#mv ${TMPDIR}/${sample}_dd.bam $dir4
