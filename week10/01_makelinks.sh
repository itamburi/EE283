
#### Make symlinks to access LRLR-KO RNAseq fastas ####

SourceDir="/dfs6/pub/itamburi/usftp21.novogene.com/01.RawData"
DestDir="/dfs6/pub/itamburi/ee283/week10"

mkdir $DestDir/datalinks

# get directory names and remove tailing slashes
ls -d "$SourceDir"/n*/ | sed 's:/*$::' > ${DestDir}/novo.ids.txt

File="${DestDir}/novo.ids.txt"

while read p
do
   echo "${p}"
   Sample=$(basename $p)
   READ1=$(find ${p}/ -type f -iname "${Sample}_1.fq.gz")
   READ2=$(find ${p}/ -type f -iname "${Sample}_2.fq.gz")
      
   for f in $READ1 $READ2
        do
           ln -s $f $DestDir/datalinks/$(basename $f)
        done

done < $File




