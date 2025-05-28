# EE283

#### Note from 5/28/25
Repository was origionally hosted locally on hpc3.rcic.uci.edu at dir `/dfs6/pub/itamburi/ee283`.
Raw data lived on a shared volume for the class, and was symlinked to to run our pipelines.
We still generated larger intermediary BAM, bigwig, etc. files. If they exceeded 100MB they were added to the .gitignore with the comand below.
On 5/28 I had to clean up my pub dir to download other data and deleted this local repo. Intermediary files can be reproduced in the future since all the code is tracable and reproducable


#### command to .gitignore files larger than 100mb upload limit
`find . -size +100M -exec echo {} \; | sed 's|^\./||' >>  .gitignore`



