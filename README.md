# ProGenFixer

ProGenFixer is a tool for detecting and correcting variants in prokaryotic genomes without the need for read mapping. It uses k-mer based analysis to identify variations between a reference genome and NGS data, and can automatically fix the reference genome based on the detected variants.

## Features

- Fast variant detection without read mapping
- Supports substitutions, insertions, and deletions
- Automatic reference genome correction
- Multi-threaded k-mer counting for improved performance
- Works with various input formats (FASTA, FASTQ, gzipped files)

## Installation

### Prerequisites

- GCC compiler (version 5.0 or higher recommended)
- zlib development libraries
- Make

### Compilation

1. Clone the repository:
   ```
   git clone https://github.com/username/ProGenFixer.git
   cd ProGenFixer
   ```

2. Compile the program:
   ```
   make
   ```

This will create the `ProGenFixer` executable in the current directory.

## Usage

### Basic Usage

```
ProGenFixer [options] Reference NGS_files > output.vcf
```

Where:
- `Reference` is the path to your reference genome in FASTA format
- `NGS_files` are one or more NGS data files (FASTQ/FASTA, can be gzipped)
- The output is directed to a VCF file

### Options

- `-k INT` : k-mer size (default: 31)
- `-c INT` : minimal k-mer coverage for variant calling (default: 3)
- `-l INT` : maximal assembly length (default: 1000)
- `-t INT` : number of threads for k-mer counting (default: 3)
- `--fix [FILE]` : Correct reference genome using detected variants (default output: fixed_reference.fna)

### Examples

1. Basic variant detection:
   ```
   ProGenFixer ref_genome.fa reads.fq > variants.vcf
   ```

2. Using paired-end reads:
   ```
   ProGenFixer ref_genome.fa reads1.fq reads2.fq > variants.vcf
   ```

3. Using multiple read files:
   ```
   ProGenFixer ref_genome.fa reads1.fq reads2.fq reads3.fq > variants.vcf
   ```

4. Fixing the reference genome:
   ```
   ProGenFixer --fix ref_genome.fa reads.fq > variants.vcf
   ```

5. Fixing the reference genome with a custom output file:
   ```
   ProGenFixer --fix=corrected_genome.fa ref_genome.fa reads.fq > variants.vcf
   ```

## Output Format

ProGenFixer outputs variants in VCF format. The columns are:

1. CHROM: Chromosome/contig name
2. POS: Position in the reference (1-based)
3. ID: Variant identifier (usually ".")
4. REF: Reference allele
5. ALT: Alternative allele
6. k-mer coverage: Coverage of the variant k-mer
7. Type: Type of variant (SUB, INS, DEL)
8. Additional info: Extra information about the variant

## How It Works

ProGenFixer uses a k-mer based approach to identify variants:

1. It builds a k-mer database from the NGS reads
2. It analyzes the reference genome to find regions where k-mer coverage drops
3. It identifies the specific variants by comparing k-mer patterns
4. It can apply these variants to create a corrected reference genome

This approach is particularly effective for prokaryotic genomes, which are typically haploid and have less complex variation patterns than eukaryotic genomes.

## Limitations

- Designed primarily for haploid genomes
- May not detect complex structural variations
- Performance depends on the chosen k-mer size and coverage thresholds

## Citation

If you use ProGenFixer in your research, please cite:

Song, L. (2025). ProGenFixer: A k-mer based tool for detecting and correcting variants in prokaryotic genomes.

## License

[License information]

## Contact

For questions or support, please contact Lifu Song at songlf@tib.cas.cn