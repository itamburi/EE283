#!/bin/sh
#SBATCH -A mseldin_lab
#SBATCH --job-name=STAR    # Job name                  
#SBATCH --cpus-per-task 16
#SBATCH --array=1
#SBATCH --mem-per-cpu=8G
#SBATCH --constraint=nvme
#SBATCH --output=/dfs6/pub/itamburi/ee283/week10/%x_%j.out   # log based on jobname
#SBATCH --error=/dfs6/pub/itamburi/ee283/week10/%x_%j.err    # error based on jobname


# source my conda software location
# source /data/homezvol0/itamburi/.bashrc
source /opt/apps/miniconda3/24.9.2/etc/profile.d/conda.sh

# activate conda enviornment with htseq
conda activate htseq

# load modules
module load star/2.7.10a

# set up enviornment variables
dir="/dfs6/pub/itamburi/ee283/week10"
prefix=`cat ${dir}/prefixes.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
#reads=$(printf '%s ' ${dir}/datalinks/${prefix}_[12])

r1="/dfs6/pub/itamburi/usftp21.novogene.com/01.RawData/n1/n1_1.fq.gz"
r2="/dfs6/pub/itamburi/usftp21.novogene.com/01.RawData/n1/n1_2.fq.gz"

cd ${TMPDIR}

STAR --genomeDir ${dir}/genomeidx \
	--readFilesIn ${r1} ${r2} \
        --outSAMtype BAM Unsorted \
        --readFilesCommand zcat \
        --outFileNamePrefix ${prefix}_ \
        --outTmpDir ${TMPDIR}/${SLURM_ARRAY_TASK_ID} \
        --runThreadN 10



dir2="${dir}/starout"

if [ ! -d "$dir2" ]; then
    mkdir -p $dir2
    echo "Directory created: $dir2"
else
    echo "Directory already exists: $dir2"
fi

mv ${prefix}_Log.out ${prefix}_Log.final.out ${prefix}_Log.progress.out $dir2
mv ${sample}_Aligned.out.bam $dir2
