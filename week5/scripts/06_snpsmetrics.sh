#!/bin/bash

# Define output file
output_file="/dfs6/pub/itamburi/ee283/week5/vcf_filt_metrics.txt"

# Write header to output file
echo -e "Step\tSNP_Count\tMean_DP\tMean_QUAL\tTiTv_Ratio" > $output_file

# List of VCF files (edit these to match your file names)
dir="/dfs6/pub/itamburi/ee283/week5/vcf"
vcf_files=("${dir}/all_variants.vcf.gz" "${dir}/filt1_variants.vcf.gz" "${dir}/filt2_variants.vcf.gz" "${dir}/filt3_variants.vcf.gz")
steps=("Initial VCF" "After Filter 1" "After Filter 2" "After Filter 3")

# Loop through each VCF file
for i in "${!vcf_files[@]}"; do
    vcf="${vcf_files[$i]}"
    step="${steps[$i]}"

    # Count total SNPs
    snp_count=$(bcftools view -H "$vcf" | wc -l)

    # Compute mean depth (DP)
    mean_dp=$(bcftools query -f '%DP\n' "$vcf" | awk '{sum+=$1} END {print sum/NR}')

    # Compute mean quality score (QUAL)
    mean_qual=$(bcftools query -f '%QUAL\n' "$vcf" | awk '{sum+=$1} END {print sum/NR}')

    # Compute Ti/Tv ratio
    titv_ratio=$(bcftools stats "$vcf" | grep "^TSTV" | tr "\n" " " | awk '{print $NF}')

    # Append results to the output file
    echo -e "$step\t$snp_count\t$mean_dp\t$mean_qual\t$titv_ratio" >> $output_file
done
