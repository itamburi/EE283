######  support script   ###########
## fasta_generate_regions.py
# example use:
#		python3 fasta_generate_regions.py $ref 10000000 >my_regions_10Mb.txt
#
# ADL modified an old script from freebayes (Feb 2021)
# intervals have to start at 1, and be non-overlapping for GTAK
# since VCF uses 1-based indexing....
#!/usr/bin/python
import sys
if len(sys.argv) == 1:
    print ("usage: ", sys.argv[0], " <fasta file or index file> <region size>")
    print ("generates a list of freebayes/bamtools region specifiers on stdout")
    print ("intended for use in creating cluster jobs")
    exit(1)
fasta_index_file = sys.argv[1]
if not fasta_index_file.endswith(".fai"):
    fasta_index_file = fasta_index_file + ".fai"
fasta_index_file = open(fasta_index_file)
region_size = int(sys.argv[2])
for line in fasta_index_file:
    fields = line.strip().split("\t")
    chrom_name = fields[0]
    chrom_length = int(fields[1])
    region_start = 1 
    while region_start < chrom_length:
        start = region_start
        end = region_start + region_size - 1
        if end > chrom_length:
            end = chrom_length
        print (chrom_name + ":" + str(region_start) + "-" + str(end))
        region_start = end + 1
