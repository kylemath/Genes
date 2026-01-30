#!/usr/bin/env python3
"""Create full chromosome 21 twin sequences"""

import subprocess
import sys

print("🧬 Creating FULL chromosome 21 twin sequences...")
print("📊 Size: ~48 million bases")
print("⏱️  This will take 5-10 minutes...")
print("")

# Get full chromosome reference
print("[1/4] Loading reference genome...")
ref = subprocess.check_output(
    "samtools faidx genome_data/reference/chr21.fa chr21", 
    shell=True
).decode()

ref_lines = ref.split('\n')[1:]  # Skip header
ref_seq = list(''.join(ref_lines))
print(f"✅ Reference loaded: {len(ref_seq):,} bases")

# Get ALL variants for chromosome 21
print("\n[2/4] Loading variants (this takes a moment)...")
vcf = "genome_data/vcf/chr21_all.vcf.gz"
cmd = f"bcftools query -s HG00096,HG00097 -f '%POS\\t%REF\\t%ALT\\t[%GT\\t]\\n' {vcf}"

variants_output = subprocess.check_output(cmd, shell=True).decode().strip()
variants = [line.split('\t') for line in variants_output.split('\n') if line]
print(f"✅ Loaded {len(variants):,} variants")

# Apply variants
print("\n[3/4] Applying variants to create twin sequences...")
seq1 = ref_seq.copy()
seq2 = ref_seq.copy()

applied = 0
progress_interval = len(variants) // 20  # Show progress 20 times

for idx, var in enumerate(variants):
    if idx % progress_interval == 0:
        pct = (idx / len(variants)) * 100
        print(f"  Progress: {pct:.0f}% ({idx:,}/{len(variants):,} variants)")
    
    if len(var) < 5:
        continue
    
    pos, ref_allele, alt_alleles, gt1, gt2 = var[0], var[1], var[2], var[3], var[4]
    offset = int(pos) - 1  # Convert to 0-indexed
    
    if offset < 0 or offset >= len(seq1):
        continue
    
    alts = alt_alleles.split(',')
    
    # Apply SNPs only (skip indels for simplicity)
    if len(ref_allele) == 1 and len(alts) > 0 and len(alts[0]) == 1:
        if gt1 in ['0/1', '1/0', '1|0', '0|1', '1/1', '1|1']:
            seq1[offset] = alts[0]
            applied += 1
        if gt2 in ['0/1', '1/0', '1|0', '0|1', '1/1', '1|1']:
            seq2[offset] = alts[0]
            applied += 1

print(f"✅ Applied {applied:,} SNP variants")

# Write output
print("\n[4/4] Writing sequences to file...")
output_file = "genome_data/fasta/CHR21_FULL_twins_comparison.fasta"

with open(output_file, "w") as f:
    f.write(">CHR21_FULL_HG00096 | British_Twin_1 | Complete_Chromosome_21\n")
    seq1_str = ''.join(seq1)
    for i in range(0, len(seq1_str), 60):
        f.write(seq1_str[i:i+60] + "\n")
    
    f.write("\n>CHR21_FULL_HG00097 | British_Twin_2 | Complete_Chromosome_21\n")
    seq2_str = ''.join(seq2)
    for i in range(0, len(seq2_str), 60):
        f.write(seq2_str[i:i+60] + "\n")

# Calculate stats
diffs = sum(1 for i in range(len(seq1)) if seq1[i] != seq2[i])
similarity = ((len(seq1)-diffs)/len(seq1)*100)

print(f"\n{'='*60}")
print(f"🎯 FULL CHROMOSOME 21 RESULTS")
print(f"{'='*60}")
print(f"Sequence length: {len(seq1):,} bases")
print(f"Differences: {diffs:,} positions")
print(f"Similarity: {similarity:.4f}%")
print(f"File size: {subprocess.check_output(f'du -h {output_file}', shell=True).decode().split()[0]}")
print(f"\n✨ Complete! File ready:")
print(f"   {output_file}")
print(f"\n⚠️  Note: This is a LARGE file. Your browser may take a moment to load it.")
