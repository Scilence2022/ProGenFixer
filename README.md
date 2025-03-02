# ProGenFixer

ProGenFixer is a tool for detecting and correcting variants in prokaryotic genomes without the need for read mapping. It uses k-mer based analysis to identify variations between a reference genome and NGS data, and can automatically fix the reference genome based on the detected variants.

## Features

- Fast variant detection without read mapping
- Supports substitutions, insertions, and deletions
- Automatic reference genome correction with iterative refinement
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
   git clone https://github.com/Scilence2022/ProGenFixer.git
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
ProGenFixer [options] Reference NGS_files
```

Where:
- `Reference` is the path to your reference genome in FASTA format
- `NGS_files` are one or more NGS data files (FASTQ/FASTA, can be gzipped)
- Output VCF files are specified by the `-o` option

### Options

- `-k INT` : k-mer size (default: 31)
- `-c INT` : minimal k-mer coverage for variant calling (default: 3)
- `-l INT` : maximal assembly length (default: 1000)
- `-t INT` : number of threads for k-mer counting (default: 3)
- `-o STR` : base name for output files (required)
- `-n INT` : number of correction iterations (default: 2)
- `--fix [FILE]` : Correct reference genome using detected variants (default output: fixed_reference.fna)

### Examples

1. Basic variant detection:
   ```
   ProGenFixer -o output ref_genome.fa reads.fq
   ```

2. Using paired-end reads:
   ```
   ProGenFixer -o output ref_genome.fa reads1.fq reads2.fq
   ```

3. Using multiple read files:
   ```
   ProGenFixer -o output ref_genome.fa reads1.fq reads2.fq reads3.fq
   ```

4. Fixing the reference genome:
   ```
   ProGenFixer -o output --fix ref_genome.fa reads.fq
   ```

5. Fixing the reference genome with a custom output file:
   ```
   ProGenFixer -o output --fix=corrected_genome.fa ref_genome.fa reads.fq
   ```

6. Using multiple correction iterations:
   ```
   ProGenFixer -o output -n 3 --fix ref_genome.fa reads.fq
   ```

## Output Files

ProGenFixer generates multiple output files:

1. VCF files for each iteration: `<output_base>.iter<N>.vcf`
2. Corrected reference files for each iteration: `<output_base>.iter<N>.fasta`
3. Final corrected reference genome after all iterations

Where `<output_base>` is the value specified with `-o` option, and `<N>` is the iteration number.

## VCF Format

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
4. It applies these variants to create a corrected reference genome
5. It can perform multiple iterations of correction to refine the results

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