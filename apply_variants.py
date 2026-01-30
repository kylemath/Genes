#!/usr/bin/env python3
"""Apply variants to create actual twin sequences"""

import subprocess

# Get reference
ref = subprocess.check_output("samtools faidx genome_data/reference/chr17.fa chr17:43044295-43049295", shell=True).decode()
ref_seq = list(''.join(ref.split('\n')[1:]))

# Get variants with genotypes
vcf = "genome_data/vcf/BRCA1_region.vcf.gz"
cmd = f"bcftools query -r 17:43044295-43049295 -s HG00096,HG00097 -f '%POS\\t%REF\\t%ALT\\t[%GT\\t]\\n' {vcf}"
variants = [line.split('\t') for line in subprocess.check_output(cmd, shell=True).decode().strip().split('\n') if line]

print(f"Processing {len(variants)} variants...")

# Create two sequences
seq1 = ref_seq.copy()
seq2 = ref_seq.copy()

applied = 0
for var in variants:
    if len(var) < 5:
        continue
    pos, ref_allele, alt_alleles, gt1, gt2 = var[0], var[1], var[2], var[3], var[4]
    
    # Convert position to 0-indexed offset
    offset = int(pos) - 43044295
    if offset < 0 or offset >= len(seq1):
        continue
    
    alts = alt_alleles.split(',')
    
    # Apply to twin 1
    if gt1 in ['0/1', '1/0', '1|0', '0|1'] and len(alts) > 0 and len(alts[0]) == 1 and len(ref_allele) == 1:
        seq1[offset] = alts[0]
        applied += 1
    elif gt1 in ['1/1', '1|1'] and len(alts) > 0 and len(alts[0]) == 1 and len(ref_allele) == 1:
        seq1[offset] = alts[0]
        applied += 1
    
    # Apply to twin 2
    if gt2 in ['0/1', '1/0', '1|0', '0|1'] and len(alts) > 0 and len(alts[0]) == 1 and len(ref_allele) == 1:
        seq2[offset] = alts[0]
        applied += 1
    elif gt2 in ['1/1', '1|1'] and len(alts) > 0 and len(alts[0]) == 1 and len(ref_allele) == 1:
        seq2[offset] = alts[0]
        applied += 1

print(f"✅ Applied {applied} simple SNP variants")

# Write output
with open("genome_data/fasta/BRCA1_twins_comparison.fasta", "w") as f:
    f.write(">BRCA1_HG00096 | British_Twin_1 | 5kb_region\n")
    seq1_str = ''.join(seq1)
    for i in range(0, len(seq1_str), 60):
        f.write(seq1_str[i:i+60] + "\n")
    
    f.write("\n>BRCA1_HG00097 | British_Twin_2 | 5kb_region\n")
    seq2_str = ''.join(seq2)
    for i in range(0, len(seq2_str), 60):
        f.write(seq2_str[i:i+60] + "\n")

# Count differences
diffs = sum(1 for i in range(len(seq1)) if seq1[i] != seq2[i])
print(f"🧬 Sequences differ at {diffs} positions")
print(f"📊 Similarity: {((len(seq1)-diffs)/len(seq1)*100):.4f}%")
print(f"\n✨ Done! Load this file in your app:")
print(f"   genome_data/fasta/BRCA1_twins_comparison.fasta")
