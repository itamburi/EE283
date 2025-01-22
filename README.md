# EE283


# command to .gitignore files larger than 100mb upload limit
`find . -size +100M -exec echo {} \; | sed 's|^\./||' >>  .gitignore`
