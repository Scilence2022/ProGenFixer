#!/bin/bash

# Check input parameters
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <ref_genome.fasta> <genome2.fasta> [output_vcf_file.vcf]"
    exit 1
fi

# Input files
ref_genome=$1
genome2=$2

# Output VCF file name (default is output_variants.vcf)
output_vcf=${3:-output_variants.vcf}

# Check if input files exist
if [ ! -f "$ref_genome" ]; then
    echo "Error: Reference genome file '$ref_genome' does not exist!"
    exit 1
fi

if [ ! -f "$genome2" ]; then
    echo "Error: Genome file '$genome2' does not exist!"
    exit 1
fi

# Create temporary directory
temp_dir=$(mktemp -d)
echo "Temporary directory created at $temp_dir"

minimap2 -a "$ref_genome" "$genome2" > "$temp_dir/alignment.sam"

# 1. Index the reference genome
#echo "Indexing reference genome..."
#bwa index "$ref_genome"

# 2. Align genome2 to the reference genome using BWA
#echo "Aligning genome2 to the reference genome using BWA..."
#bwa mem "$ref_genome" "$genome2" > "$temp_dir/alignment.sam"
# 3. Convert SAM to BAM format and sort
echo "Converting SAM to BAM and sorting..."
samtools view -Sb "$temp_dir/alignment.sam" | samtools sort -o "$temp_dir/alignment.sorted.bam"

# 4. Index the sorted BAM file
echo "Indexing sorted BAM file..."
samtools index "$temp_dir/alignment.sorted.bam"

# 5. Call variants using bcftools
echo "Calling variants using bcftools..."
bcftools mpileup -f "$ref_genome" "$temp_dir/alignment.sorted.bam" | bcftools call -mv --ploidy 1 -Ob -o "$temp_dir/variants.bcf"

# 6. Convert BCF to VCF format
echo "Converting BCF to VCF..."
bcftools view "$temp_dir/variants.bcf" > "$temp_dir/variants.vcf"

# 7. Move the final VCF file to the target path
echo "Saving VCF file to '$output_vcf'..."
mv "$temp_dir/variants.vcf" "$output_vcf"

# Clean up temporary files
rm -rf "$temp_dir"

echo "Comparison complete. Variant report saved to '$output_vcf'"
