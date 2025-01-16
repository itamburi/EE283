#class data directory: /data/class/ecoevo283/public/Bioinformatics_Course

#### DNA-SEQ
SourceDir="/data/class/ecoevo283/public/Bioinformatics_Course/DNAseq"
DestDir="/dfs6/pub/itamburi/ee283/datalinks/DNAseq"
FILES="$SourceDir/*"
for f in $FILES
do
   ff=$(basename $f)
   echo "Processing $ff file..."
   ln -s $SourceDir/$ff $DestDir/$ff
done



#### ATAC-SEQ
SourceDir="/data/class/ecoevo283/public/Bioinformatics_Course/ATACseq"
DestDir="/dfs6/pub/itamburi/ee283/datalinks/ATACseq"

tail -n +2 ${SourceDir}/README.ATACseq.txt | head -n -3 > ${DestDir}/ATACseq.labels.txt
File="${DestDir}/ATACseq.labels.txt"

# look at the end of the loop, p in this context assigns one line at a time to be read
while read p
do
   echo "${p}"
   barcode=$(echo $p | cut -f1 -d" ")
   genotype=$(echo $p | cut -f2 -d" ")
   tissue=$(echo $p | cut -f3 -d" ")
   bioRep=$(echo $p | cut -f4 -d" ")
   READ1=$(find ${SourceDir}/ -type f -iname "*_${barcode}_R1.fq.gz")
   READ2=$(find ${SourceDir}/ -type f -iname "*_${barcode}_R2.fq.gz")
   
   sample="${barcode}_${genotype}_${tissue}_${bioRep}"
   
   for f in $READ1 $READ2
	do
   	   ff=$(basename $f)
   	   echo "Processing $ff file..."
   	   read=$(echo "$ff" | grep -o 'R[12]')
	   ln -s $f $DestDir/${sample}_$read
	done

done < $File



#### RNA-SEQ
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

# final format of the link will be sampleID (from col 12 "FullSampleName" in metadata) followed by _read
       done

done < $File

