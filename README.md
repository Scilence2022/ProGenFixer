# ProGenFixer

ProGenFixer is a super-fast and accurate tool for detecting and correcting variants in prokaryotic genomes without the need for read mapping. It uses k-mer based analysis to identify variations between a reference genome and NGS data, and can automatically fix the reference genome based on the detected variants.

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
- `-c INT` : minimal k-mer coverage for detection of variation regions (default: 2)
- `-a INT` : minimal k-mer coverage for variant calling (default: 5)
- `-l INT` : maximal assembly length for variant calling (default: 1000)
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

ProGenFixer outputs variants in VCF format (VCFv4.5 compatible). The VCF files include the following header information:

```
##fileformat=VCFv4.5
##ProGenFixerVersion=v1.0
##ProGenFixerCommand=[command line used]
##INFO=<ID=KMER_COV,Number=1,Type=Integer,Description="K-mer coverage of the variant path">
##INFO=<ID=VARTYPE,Number=1,Type=String,Description="Variant type (INS, DEL, SUB)">
```

The VCF records contain the standard columns:

1. CHROM: Chromosome/contig name
2. POS: Position in the reference (1-based)
3. ID: Variant identifier (always "." in current implementation)
4. REF: Reference allele
5. ALT: Alternative allele
6. QUAL: Quality score (always "." in current implementation)
7. FILTER: Filter status (always "PASS" in current implementation)
8. INFO: Additional information in the format `KMER_COV=<value>;VARTYPE=<type>`

Where `<type>` is one of:
- SUB: Substitution
- INS: Insertion
- DEL: Deletion

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

Song et al. (2025). ProGenFixer: an ultra-fast and accurate tool for correcting prokaryotic genome sequences using a mapping-free algorithm

## License

ProGenFixer is licensed under the GNU General Public License v3.0 (GPL-3.0).

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

## Contact

For questions or support, please contact Lifu Song at songlf@tib.cas.cn
