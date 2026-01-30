#!/usr/bin/env python3
"""Create twin sequences from chromosome 21"""

import subprocess
import sys

print("🧬 Creating chromosome 21 twin sequences...")
print("📊 This is ~48 million bases (16,000x larger than BRCA1 sample!)")
print("")

# Extract first 100kb region (manageable for analysis)
region_start = 5000000
region_size = 100000
region_end = region_start + region_size

print(f"Extracting region: chr21:{region_start}-{region_end} ({region_size:,} bases)")

# Get reference
ref = subprocess.check_output(
    f"samtools faidx genome_data/reference/chr21.fa chr21:{region_start}-{region_end}", 
    shell=True
).decode()
ref_seq = list(''.join(ref.split('\n')[1:]))

print(f"✅ Reference extracted: {len(ref_seq):,} bases")

# Get variants
vcf = "genome_data/vcf/chr21_all.vcf.gz"
cmd = f"bcftools query -r 21:{region_start}-{region_end} -s HG00096,HG00097 -f '%POS\\t%REF\\t%ALT\\t[%GT\\t]\\n' {vcf}"
variants = [line.split('\t') for line in subprocess.check_output(cmd, shell=True).decode().strip().split('\n') if line]

print(f"✅ Found {len(variants):,} variants in region")

# Create sequences
seq1 = ref_seq.copy()
seq2 = ref_seq.copy()

applied = 0
for var in variants:
    if len(var) < 5:
        continue
    pos, ref_allele, alt_alleles, gt1, gt2 = var[0], var[1], var[2], var[3], var[4]
    
    offset = int(pos) - region_start
    if offset < 0 or offset >= len(seq1):
        continue
    
    alts = alt_alleles.split(',')
    
    # Apply SNPs only (for simplicity)
    if len(ref_allele) == 1 and len(alts) > 0 and len(alts[0]) == 1:
        if gt1 in ['0/1', '1/0', '1|0', '0|1', '1/1', '1|1']:
            seq1[offset] = alts[0]
            applied += 1
        if gt2 in ['0/1', '1/0', '1|0', '0|1', '1/1', '1|1']:
            seq2[offset] = alts[0]
            applied += 1

print(f"✅ Applied {applied:,} SNP variants")

# Write output
with open("genome_data/fasta/CHR21_twins_comparison.fasta", "w") as f:
    f.write(">CHR21_HG00096 | British_Twin_1 | 100kb_region\n")
    seq1_str = ''.join(seq1)
    for i in range(0, len(seq1_str), 60):
        f.write(seq1_str[i:i+60] + "\n")
    
    f.write("\n>CHR21_HG00097 | British_Twin_2 | 100kb_region\n")
    seq2_str = ''.join(seq2)
    for i in range(0, len(seq2_str), 60):
        f.write(seq2_str[i:i+60] + "\n")

diffs = sum(1 for i in range(len(seq1)) if seq1[i] != seq2[i])
print(f"\n🎯 Results:")
print(f"   Sequence length: {len(seq1):,} bases")
print(f"   Differences: {diffs:,} positions")
print(f"   Similarity: {((len(seq1)-diffs)/len(seq1)*100):.4f}%")
print(f"\n✨ Done! File ready:")
print(f"   genome_data/fasta/CHR21_twins_comparison.fasta")
print(f"\n📈 Comparison:")
print(f"   BRCA1 sample: 5,000 bases")
print(f"   This sample: {len(seq1):,} bases (20x larger!)")
print(f"   Full chr21: ~48,000,000 bases (480x larger than this!)")
