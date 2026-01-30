# 🧬 Intelligent Variant Classification Guide

## 🎯 Smart Grouping of Adjacent Differences

The app now intelligently groups adjacent sequence differences into meaningful genetic variants, just like professional bioinformatics tools!

---

## 🔬 Variant Types

### **SNP - Single Nucleotide Polymorphism**
```
Position: 125
Twin A: A
Twin B: G

Classification: SNP (most common)
Impact: May be synonymous or missense
```
**Example:** One base differs

### **DNP - Di-Nucleotide Polymorphism**
```
Position: 250-251
Twin A: AT
Twin B: GC

Classification: DNP (rare)
Impact: Affects 2 consecutive bases
```
**Example:** Two adjacent bases differ

### **TNP - Tri-Nucleotide Polymorphism**
```
Position: 380-382
Twin A: ATG
Twin B: GCA

Classification: TNP (one codon)
Impact: Complete codon change
```
**Example:** Three adjacent bases differ (whole codon)

### **MNP - Multi-Nucleotide Polymorphism**
```
Position: 500-505
Twin A: AAGATG
Twin B: GCTTAA

Classification: MNP (6 bases)
Impact: Multiple codons affected
```
**Example:** 4+ adjacent bases differ

---

## 🧠 Why This Matters

### Before (Naive Approach):
```
❌ Difference 1: Position 500 (A→G)
❌ Difference 2: Position 501 (A→C)
❌ Difference 3: Position 502 (G→T)
❌ Difference 4: Position 503 (A→T)
❌ Difference 5: Position 504 (T→A)
❌ Difference 6: Position 505 (G→A)

Result: Looks like 6 separate mutations! 😵
```

### After (Smart Grouping):
```
✅ Variant 1: MNP at position 500-505
   AAGATG → GCTTAA
   
Result: One complex mutation event! 🎯
```

---

## 📊 What You'll See Now

### Variant Summary Card
```
🔍 Genetic Variants Detected
Found 3 genetic variants: 2 SNPs, 1 MNP

┌─────────┬─────────┬─────────┐
│  2 SNP  │  0 DNP  │  1 MNP  │
└─────────┴─────────┴─────────┘
```

### Individual Variant Cards

#### SNP Example:
```
┌─────────────────────────────────────────┐
│ [SNP] Variant 1                         │
│ Position: 125                           │
│                                         │
│ Single Nucleotide Polymorphism          │
│                                         │
│    A  →  G                              │
│                                         │
│ Context: ...ATCGA[A]TCGAT...           │
│          ...ATCGA[G]TCGAT...           │
│                                         │
│ 📍 Affects codon #42 (frame 2)         │
│ ⚠️ In-frame (maintains reading frame)   │
└─────────────────────────────────────────┘
```

#### MNP Example:
```
┌─────────────────────────────────────────┐
│ [MNP] Variant 2                         │
│ Position: 500-505                       │
│                                         │
│ Multi-Nucleotide Polymorphism (6 bases)│
│                                         │
│  A A G A T G  →  G C T T A A          │
│                                         │
│ Context: ...CTGAA[AAGATG]CCTATGC...   │
│          ...CTGAA[GCTTAA]CCTATGC...   │
│                                         │
│ 📍 Affects codon #167 (frame 1)        │
│ ⚠️ In-frame (maintains reading frame)   │
└─────────────────────────────────────────┘
```

---

## 🎨 Visual Classification

### Color-Coded Badges:
- **🟡 SNP** - Yellow/Orange (warning-color)
- **🟠 DNP** - Orange (#f97316)
- **🔴 TNP** - Red (#ef4444)
- **🔴 MNP** - Dark Red (#dc2626)

### Border Colors Match Badges:
- Left border of variant card matches the variant type
- Darker colors = more bases affected
- Instantly see severity!

---

## 💡 Biological Significance

### SNPs (Single Base)
**Clinical Relevance:**
- Most common type of genetic variation
- Can affect protein function
- Used in GWAS studies
- Pharmacogenomics markers

**Example:** BRCA1 185delAG

### TNPs (Three Bases = One Codon)
**Clinical Relevance:**
- Affects exactly one amino acid
- Common in protein evolution
- Conservative vs. non-conservative changes matter

**Example:** Sickle cell (GAG→GTG in HBB)

### MNPs (Multiple Bases)
**Clinical Relevance:**
- Complex mutation events
- May indicate:
  - Gene conversion
  - Recombination
  - Multiple simultaneous mutations
  - Sequencing artifacts (check carefully!)

**Example:** Deletions in cystic fibrosis

---

## 🔍 Detailed Variant Information

### What Each Field Means:

#### **Type & Description**
```
MNP - Multi-Nucleotide Polymorphism (6 bases)
```
- Type: Classification
- Description: Human-readable explanation

#### **Position Range**
```
Position: 500-505
```
- Start position (1-indexed)
- End position (for multi-base variants)

#### **Base Changes**
```
AAGATG → GCTTAA
```
- Original sequence (Twin A)
- Changed sequence (Twin B)
- Each base color-coded

#### **Context Display**
```
...CTGAA[AAGATG]CCTATGC...
```
- 10 bases before and after
- Variant highlighted in brackets
- Full context for understanding location

#### **Codon Information**
```
📍 Affects codon #167 (frame 1)
```
- Which codon number is affected
- Reading frame (1, 2, or 3)
- Useful for predicting protein changes

#### **Impact Assessment**
```
⚠️ In-frame (maintains reading frame)
```
- **In-frame**: Length divisible by 3, reading frame maintained
- **Frameshift possible**: Length not divisible by 3, may shift reading

---

## 🧪 Real-World Examples

### Example 1: BRCA1 Pathogenic SNP
```
Type: SNP
Position: 5382
Change: C → T
Codon: 1755 (frame 3)
Impact: Nonsense mutation (creates STOP)
Clinical: Highly pathogenic for breast/ovarian cancer
```

### Example 2: CFTR Common Deletion
```
Type: MNP (actually a 3-base deletion)
Position: 1521-1523
Change: CTT → ---
Codon: 508
Impact: Frameshift, ΔF508
Clinical: Most common CF mutation
```

### Example 3: Hemoglobin S (Sickle Cell)
```
Type: SNP
Position: 17
Change: A → T
Codon: 6
Impact: Glu → Val (GAG → GTG)
Clinical: Causes sickle cell disease when homozygous
```

---

## 📈 Statistics You'll See

### Variant Distribution
```
Total Variants: 5
├── SNPs: 3 (60%)
├── DNPs: 1 (20%)
├── TNPs: 0 (0%)
└── MNPs: 1 (20%)
```

### Sequence Identity
```
99.8718% identical
5 variants across 780 bases
1 variant per 156 bases (average)
```

---

## 🎯 How Grouping Works (Algorithm)

### Step 1: Find All Differences
```javascript
// Compare position by position
seq1: ATGAAGATGC
seq2: ATGGCTTAAC
      ^^^    ^    
Diffs at: 3,4,5,9
```

### Step 2: Group Adjacent Positions
```javascript
// Check if positions are consecutive
Position 3,4,5 → Group as MNP (3 bases)
Position 9     → Separate SNP (1 base)
```

### Step 3: Classify by Length
```javascript
if (length === 1) → SNP
if (length === 2) → DNP
if (length === 3) → TNP
if (length >= 4) → MNP
```

### Step 4: Analyze Impact
```javascript
// Check reading frame
if (length % 3 === 0) → In-frame
else → Possible frameshift
```

---

## 🔬 Advanced Features

### Frame Analysis
Every variant shows which codon and reading frame it affects:
- **Codon number**: Position ÷ 3 + 1
- **Frame**: Position mod 3 (1, 2, or 3)

### Context Windows
- Shows ±10 bases around variant
- Helps understand local sequence context
- Important for primer design

### Impact Prediction
- **In-frame**: Maintains codon structure
- **Frameshift**: May disrupt downstream codons
- Basis for further protein analysis

---

## 💻 Technical Implementation

### Grouping Algorithm (Pseudocode)
```python
def group_variants(differences):
    variants = []
    i = 0
    
    while i < len(differences):
        start = differences[i]
        end = start
        
        # Extend while consecutive
        while (i+1 < len(differences) and 
               differences[i+1] == differences[i] + 1):
            i += 1
            end = differences[i]
        
        # Create variant
        length = end - start + 1
        type = classify(length)
        variants.append(Variant(start, end, type))
        
        i += 1
    
    return variants
```

---

## 🎓 Learning Guide

### For Beginners:
1. Load twin sample data
2. Look at variant summary
3. Click on each variant
4. Read the descriptions
5. Understand SNPs first, then MNPs

### For Intermediate:
1. Compare different genes
2. Look for patterns (SNP clusters)
3. Analyze codon impacts
4. Predict protein changes

### For Advanced:
1. Download real clinical data
2. Compare to ClinVar database
3. Assess pathogenicity
4. Design targeted therapies

---

## 📚 Related Concepts

### Variant Nomenclature (HGVS)
Standard notation used in genetics:
```
c.500A>G        (SNP at position 500)
c.500_505delins (MNP, deletion-insertion)
p.Arg506Gly     (protein level change)
```

### VCF Format
Standard file format for variants:
```
#CHROM  POS  ID    REF  ALT  QUAL
chr17   500  rs123 A    G    99
```

### Clinical Databases
- **ClinVar**: Clinical significance
- **dbSNP**: SNP repository
- **COSMIC**: Cancer mutations
- **gnomAD**: Population frequencies

---

## 🆘 Troubleshooting

### "Too many variants detected"
**Check:**
- Are you comparing the same gene region?
- Different species?
- Sequencing quality issues?

**Expected:** 1-5 variants for twins in same gene

### "All variants show as SNPs"
**This means:** No adjacent differences found
**Is normal for:** Most real twin comparisons

### "MNP with many bases"
**Could indicate:**
- Recombination event
- Gene conversion
- Technical artifact (verify!)
- Different gene paralogs

---

## ✅ Best Practices

1. **Always check context** - Look at surrounding sequence
2. **Verify MNPs** - Large MNPs are rare, may be errors
3. **Consider reading frame** - In-frame vs frameshift matters
4. **Look up in databases** - Check if variant is known
5. **Replicate findings** - Verify with independent sequencing

---

## 🎯 Quick Reference

| Variant | Bases | Common? | Clinical Relevance |
|---------|-------|---------|-------------------|
| SNP     | 1     | ✅ Very | High - most studied |
| DNP     | 2     | ⚠️ Rare | Medium - less understood |
| TNP     | 3     | ⚠️ Rare | High - affects one codon |
| MNP     | 4+    | ❌ Very rare | High - complex events |

---

**The app now thinks like a geneticist! 🧬🧠**

**Try loading the twin sample to see the intelligent grouping in action!**

