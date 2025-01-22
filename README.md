# EE283


# command to .gitignore files larger than 100mb upload limit
`find . -size +100M -exec echo {} \; | sed 's|^\./||' >>  .gitignore`



# 01/22: prob2:
# > Check that trimmomatic was run with the correct parameters, or if the parameters could be improved
