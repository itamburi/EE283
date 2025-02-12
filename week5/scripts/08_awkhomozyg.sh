#!/bin/bash

# To me, this problem made much more sense to perform in R (see scripts/plot_biallelic.R)
# In R it is super easy to transpose, filter 1's, calculate the sums across cols, and plot
# ... However, to demosntrate this in awk, we could also perform the following:
# X_1Mb_snps.012 has N cols for bp positions, and rows for our genotypes. This is cumbersome to compute sums on.
# ... so first I opened are and did a simple transpose with t() and made the file X_1Mb_snps.012_transposed
# the following awk will work to produce the desired result that we can plot in R


dir="/dfs6/pub/itamburi/ee283/week5"
awk '
{
    has_one = 0;
    sum = 0;
    
    for (i = 1; i <= NF; i++) {
        if ($i == 1) has_one = 1;
        sum += $i;
    }
    
    if (!has_one && sum == 4) {
        print $0;
    }
}' ${dir}/X_1Mb_snps.012_transposed > ${dir}/X_1Mb_snps.012_filt
