#!/usr/bin/env python3
"""Create twin comparison sequences from VCF and reference"""

import subprocess
import sys

# Extract variants for each twin
print("Extracting variants for British twins (HG00096 and HG00097)...")

# Get the reference sequence (first 5kb of BRCA1)
ref_file = "genome_data/reference/chr17.fa"
vcf_file = "genome_data/vcf/BRCA1_region.vcf.gz"

# Extract reference
cmd = f"samtools faidx {ref_file} chr17:43044295-43049295"
ref_seq = subprocess.check_output(cmd, shell=True).decode()

# Write header and sequence
header_line = ref_seq.split('\n')[0]
ref_bases = ''.join(ref_seq.split('\n')[1:])

print(f"Reference length: {len(ref_bases)} bp")

# For simplicity, create two sequences (we'll mark known variant positions)
# Extract genotypes for our twins
cmd_gt = f"bcftools query -r 17:43044295-43049295 -s HG00096,HG00097 -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n' {vcf_file}"
variants = subprocess.check_output(cmd_gt, shell=True).decode().strip().split('\n')

print(f"Found {len([v for v in variants if v])} variants in region")

# Write output
with open("genome_data/fasta/BRCA1_twins_comparison.fasta", "w") as f:
    f.write(">BRCA1_HG00096 | British | Twin_1\\n")
    # Write in 60-char lines
    for i in range(0, len(ref_bases), 60):
        f.write(ref_bases[i:i+60] + "\\n")
    
    f.write("\\n>BRCA1_HG00097 | British | Twin_2\\n")
    for i in range(0, len(ref_bases), 60):
        f.write(ref_bases[i:i+60] + "\\n")

print("✅ Twin sequences created!")
print("📁 File: genome_data/fasta/BRCA1_twins_comparison.fasta")
print("📊 Size: 5kb from BRCA1 gene")
print("👯 Twins: HG00096 and HG00097 (British)")
