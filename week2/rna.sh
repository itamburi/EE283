# /data/class/ecoevo283/public/Bioinformatics_Course


# RNA-SEQ
SourceDir="/data/class/ecoevo283/public/Bioinformatics_Course/RNAseq"
DestDir="/dfs6/pub/itamburi/ee283/datalinks/RNAseq"

tail -n +2 ${SourceDir}/RNAseq384_SampleCoding.txt > ${DestDir}/RNAseq.labels.txt
File="${DestDir}/RNAseq.labels.txt"

# look at the end of the loop, p in this context assigns one line at a time to be read
while read p
do
#   echo "${p}"
   plex=$(echo $p | cut -f2 -d" " | sed 's/_[^_]*$//')
   plex="Project_$plex"

   number=$(echo $p | cut -f1 -d" ")
   number="Sample_$number"

   barcode=$(echo $p | cut -f4 -d" ")
   sampleid=$(echo $p | cut -f12 -d" ")
   
   path="${SourceDir}/RNAseq384plex_flowcell01/${plex}/${number}"

   READ1=$(find ${path}/ -type f -iname "*_${barcode}_*_R1_*.fastq.gz")
   READ2=$(find ${path}/ -type f -iname "*_${barcode}_*_R2_*.fastq.gz")
   
   for f in $READ1 $READ2
       do
          ff=$(basename $f)
          echo "Processing $ff file..."
          read=$(echo "$ff" | grep -o 'R[12]')
          echo $read
          ln -s $f $DestDir/${sampleid}_$read
       done

done < $File

