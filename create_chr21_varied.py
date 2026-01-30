#!/usr/bin/env python3
import subprocess

# Try APP gene region (Alzheimer's related, chr21:25880550-26171128)
region_start = 25880550
region_size = 50000  # 50kb
region_end = region_start + region_size

print(f"🧬 Extracting APP gene region: chr21:{region_start}-{region_end}")

ref = subprocess.check_output(
    f"samtools faidx genome_data/reference/chr21.fa chr21:{region_start}-{region_end}", 
    shell=True
).decode()
ref_seq = list(''.join(ref.split('\n')[1:]))

vcf = "genome_data/vcf/chr21_all.vcf.gz"
cmd = f"bcftools query -r 21:{region_start}-{region_end} -s HG00096,HG00097 -f '%POS\\t%REF\\t%ALT\\t[%GT\\t]\\n' {vcf}"
variants = [line.split('\t') for line in subprocess.check_output(cmd, shell=True).decode().strip().split('\n') if line]

print(f"✅ Found {len(variants)} variants")

seq1, seq2 = ref_seq.copy(), ref_seq.copy()
applied = 0

for var in variants:
    if len(var) < 5: continue
    pos, ref_allele, alt_alleles, gt1, gt2 = var[0], var[1], var[2], var[3], var[4]
    offset = int(pos) - region_start
    if offset < 0 or offset >= len(seq1): continue
    alts = alt_alleles.split(',')
    
    if len(ref_allele) == 1 and len(alts) > 0 and len(alts[0]) == 1:
        if gt1 in ['0/1', '1/0', '1|0', '0|1', '1/1', '1|1']:
            seq1[offset] = alts[0]
            applied += 1
        if gt2 in ['0/1', '1/0', '1|0', '0|1', '1/1', '1|1']:
            seq2[offset] = alts[0]
            applied += 1

with open("genome_data/fasta/CHR21_twins_comparison.fasta", "w") as f:
    f.write(">CHR21_HG00096 | APP_gene_region | 50kb\n")
    seq1_str = ''.join(seq1)
    for i in range(0, len(seq1_str), 60):
        f.write(seq1_str[i:i+60] + "\n")
    f.write("\n>CHR21_HG00097 | APP_gene_region | 50kb\n")
    seq2_str = ''.join(seq2)
    for i in range(0, len(seq2_str), 60):
        f.write(seq2_str[i:i+60] + "\n")

diffs = sum(1 for i in range(len(seq1)) if seq1[i] != seq2[i])
print(f"✅ Applied {applied} variants")
print(f"🎯 {diffs} differences between twins")
print(f"📊 {((len(seq1)-diffs)/len(seq1)*100):.4f}% similarity")
print(f"✨ genome_data/fasta/CHR21_twins_comparison.fasta")
