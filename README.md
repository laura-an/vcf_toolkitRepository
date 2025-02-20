# **vcf_toolkit**
`vcf_toolkit.py` is a is a Python tool designed to efficiently handle and analyze VCF (Variant Call Format) files, providing various functions such 
as position-based search, row counting, sample and chromosome exploration, among others.

## Installation
This tool requires `python3`, `bcftools` and `pandas` to work correctly.

## Usage
It can be used as a module. Also, the tool can be run directly using the following command:

`./vcf_toolkit.py <command> -input file.vcf [options]`

## Example
`./vcf_toolkit.py search -input example.vcf -position chr1:10000`

## Available Commands
- `nrows` → Returns the total number of rows.
- `nrows_nh` → Returns the number of rows excluding the header.
- `search` → Searches for a specific position (chr:pos).
- `nchrom` → Returns the number of rows per chromosome.
- `chrom` → Filters rows corresponding to a specific chromosome.
- `nsamples` → Returns the total number of samples.
- `samples` → Lists the sample names.
- `filter` → Counts occurrences of PASS and FAIL in the FILTER column.
- `description` → Displays a summary of the VCF file.
- `readVcf` → Loads the VCF file into a pandas DataFrame.

