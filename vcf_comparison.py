#!/usr/bin/env python3

import subprocess
import argparse
import os

def run_bcftools_isec(vcf1, vcf2, output_dir):
    """
    Use bcftools isec to compare two VCF files and save the results in the output directory.
    """
    command = [
        'bcftools', 'isec', vcf1, vcf2, '-p', output_dir
    ]
    subprocess.run(command, check=True)

def parse_vcf_file(file_path, filter_words=None):
    """
    Parse a VCF file and count the number of variants.
    If filter_words is provided, only count lines that contain the specified words.
    
    Args:
        file_path: Path to the VCF file
        filter_words: List of words to filter by, default is None (no filtering)
    """
    count = 0
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            if filter_words:
                if any(word in line for word in filter_words):
                    count += 1
            else:
                count += 1
    
    return count

def calculate_metrics(tp, fp, fn):
    """
    Calculate Sensitivity and Precision
    """
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    return sensitivity, precision

def bgzip_and_index(vcf_file):
    """
    If the VCF file is not compressed with bgzip, compress it using bgzip and create an index.
    """
    compressed_vcf = vcf_file + '.gz'
    
    # If the file is not compressed, compress it using bgzip
    if not os.path.exists(compressed_vcf):
        print(f"Compressing {vcf_file} using bgzip...")
        try:
            # Compress the file with bgzip
            subprocess.run(['bgzip', '-k', '-f', vcf_file], check=True)
            print(f"File compressed successfully: {compressed_vcf}")
        except subprocess.CalledProcessError as e:
            print(f"Error during bgzip compression: {e}")
            raise
    
    # Create an index file
    print(f"Creating index for {compressed_vcf}...")
    try:
        subprocess.run(['tabix', '-f', '-p','vcf', compressed_vcf], check=True)
        print(f"Index created successfully: {compressed_vcf}.tbi")
    except subprocess.CalledProcessError as e:
        print(f"Error during tabix indexing: {e}")
        raise

    return compressed_vcf

def main(vcf1, vcf2, output_dir, filter_words=None):
    # Step 1: Compress VCF files and generate index files
    print("Compressing and indexing VCF files if necessary...")
    vcf1_compressed = bgzip_and_index(vcf1)
    vcf2_compressed = bgzip_and_index(vcf2)
    
    # Step 2: Compare VCF files using bcftools isec
    print("Running bcftools isec to compare the VCF files...")
    run_bcftools_isec(vcf1_compressed, vcf2_compressed, output_dir)
    
    # Step 3: Calculate TP, FP, FN
    tp_file = f"{output_dir}/0002.vcf"
    fp_file = f"{output_dir}/0001.vcf"
    fn_file = f"{output_dir}/0000.vcf"
    
    filter_msg = f" (filtered by: {', '.join(filter_words)})" if filter_words else ""
    
    tp_count = parse_vcf_file(tp_file, filter_words)
    fp_count = parse_vcf_file(fp_file, filter_words)
    fn_count = parse_vcf_file(fn_file, filter_words)
    
    print(f"True Positives (TP){filter_msg}: {tp_count}")
    print(f"False Positives (FP){filter_msg}: {fp_count}")
    print(f"False Negatives (FN){filter_msg}: {fn_count}")
    
    # Step 4: Calculate Sensitivity and Precision
    sensitivity, precision = calculate_metrics(tp_count, fp_count, fn_count)
    
    print(f"Sensitivity{filter_msg}: {sensitivity:.4f}")
    print(f"Precision{filter_msg}: {precision:.4f}")
    
    # Step 5: Output the results
    with open(f"{output_dir}/comparison_metrics.txt", "w") as out_file:
        out_file.write(f"True Positives (TP){filter_msg}: {tp_count}\n")
        out_file.write(f"False Positives (FP){filter_msg}: {fp_count}\n")
        out_file.write(f"False Negatives (FN){filter_msg}: {fn_count}\n")
        out_file.write(f"Sensitivity{filter_msg}: {sensitivity:.4f}\n")
        out_file.write(f"Precision{filter_msg}: {precision:.4f}\n")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Compare two VCF files and calculate TP, FP, FN, Sensitivity, and Precision.")
    parser.add_argument("vcf1", help="Path to the first VCF file (ground truth)")
    parser.add_argument("vcf2", help="Path to the second VCF file (predictions)")
    parser.add_argument("output_dir", help="Directory to save the comparison output")
    parser.add_argument("--filter", nargs='+', help="Optional list of words to filter VCF records (e.g., 'INDEL')")
    
    args = parser.parse_args()

    # Check if the output directory exists, if not, create it
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Call the main function to compare VCF files and compute metrics
    main(args.vcf1, args.vcf2, args.output_dir, args.filter)
