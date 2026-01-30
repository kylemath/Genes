# 🧬 Full Sequence Alignment View Guide

## 🎨 The Big Picture - See Everything at Once!

The new **Full Sequence Alignment View** gives you a complete, side-by-side comparison of entire twin sequences with all variants highlighted in context. It's like having a professional sequence alignment tool built right into your browser!

---

## ✨ What You Get

### 1. **Side-by-Side Sequence Display**
```
Position:  1                                                      60
Twin 1:    ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAA
           |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Twin 2:    ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAA

Position:  61                                                    120
Twin 1:    ATCTTAGAGTGTCCCATCTGGTAAGTCAGCACAGAGGGCCAGCTCATTCCTGCCCTCTGC
           |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Twin 2:    ATCTTAGAGTGTCCCATCTGGTAAGTCAGCACAGAGGGCCAGCTCATTCCTGCCCTCTGC
```

### 2. **Visual Variant Highlighting**
- 🟡 **SNPs**: Yellow/orange glow
- 🟠 **DNPs**: Orange glow
- 🔴 **TNPs/MNPs**: Red glow
- ⚪ **Matches**: Subtle gray background

### 3. **Match Indicators**
- `|` = Bases match perfectly
- `×` = Bases differ (variant position)

### 4. **Interactive Features**
- ✅ **Hover** over any base to highlight both twins
- ✅ **Click** any base to see detailed position info
- ✅ **Scroll** horizontally to see full length
- ✅ **Scroll** vertically for long sequences

### 5. **Variant Density Map**
Visual bar chart showing where variants cluster along the sequence!

---

## 🎯 How to Use

### Step 1: Load Twin Data
```bash
open index.html
# Click "👯 Twin BRCA1 Comparison"
# Click "🔬 Analyze Sequence"
```

### Step 2: Scroll to Alignment Section
Look for the **"🧬 Full Sequence Alignment"** card

### Step 3: Explore!
- **Read the sequence** in 60-base blocks
- **Spot variants** by their bright colors
- **Check the density map** at the bottom
- **Hover and click** for details

---

## 🖼️ Visual Layout

```
┌─────────────────────────────────────────────────────────────┐
│ 🧬 Full Sequence Alignment                                  │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│ Legend: [Match] [SNP] [DNP] [MNP/TNP]    3 variants/780bp  │
│                                                              │
│ ┌──────────────────────────────────────────────────────┐   │
│ │  1   ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATT  │   │
│ │      |||||||||||||||||||||||||||||||||||||||||||||| │   │
│ │  1   ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATT  │   │
│ │                                                       │   │
│ │  61  AATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGGTAAGTCAG  │   │
│ │      |||||||||||||||||||||||||||||||||||||||||||||| │   │
│ │  61  AATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGGTAAGTCAG  │   │
│ │                                                       │   │
│ │  ... (continues for full sequence)                    │   │
│ └──────────────────────────────────────────────────────┘   │
│                                                              │
│ Variant Density Map:                                        │
│ ┌──────────────────────────────────────────────────────┐   │
│ │ ||  ||          ||                    ||        ||    │   │
│ └──────────────────────────────────────────────────────┘   │
│ Start ─────────── Position ─────────────────────── End     │
└─────────────────────────────────────────────────────────────┘
```

---

## 🔍 Interactive Features in Detail

### Hover Over Bases
```
Action: Mouse over any base
Result: 
  - Corresponding position in both sequences highlights
  - Easy to trace exact position across twins
  - Visual guide with scale effect
```

### Click on Bases
```
Action: Click any base
Result: Tooltip appears showing:
  ┌─────────────────────┐
  │ Position 645        │
  │ Twin 1: A           │
  │ Twin 2: G           │
  │ ⚠️ Difference       │
  │ Codon: #215         │
  │ Frame: 3            │
  └─────────────────────┘
  
Auto-hides after 3 seconds
```

### Click on Density Map
```
Action: Click any bar in density map
Result:
  - Scrolls to that section of sequence
  - Shows variants in that region
  - Quick navigation for long sequences
```

---

## 📊 Variant Density Map

### What It Shows
The density map divides your sequence into 100 segments and shows how many variants appear in each region.

```
High density region:
████ = Many variants here (hotspot!)

Low density region:
█    = Few or no variants (conserved region)
```

### Why It's Useful
- **Identify hotspots**: Where most changes occur
- **Find conserved regions**: Important functional areas
- **Navigate quickly**: Click to jump to interesting regions
- **Pattern recognition**: Are variants clustered or spread?

### Example Interpretation
```
Density Map:
[     ||||||||     |     |           ||||    ]
^     ^        ^    ^     ^           ^       ^
|     |        |    |     |           |       |
Start Hotspot  End  SNP   Mid        Cluster End

Interpretation:
- Positions 200-400: High variant density (functional region?)
- Positions 600-800: Conserved (critical for function?)
- Position 450: Isolated SNP (benign variant?)
```

---

## 🎨 Color-Coding System

### Visual Hierarchy

#### **Matches (Most Common)**
```css
Background: Light gray (barely visible)
Text: White
Purpose: Show normal, identical bases
```

#### **SNPs (Common)**
```css
Background: Orange (#f59e0b)
Glow: Orange shadow
Text: Black
Purpose: Highlight single base changes
```

#### **DNPs (Rare)**
```css
Background: Dark Orange (#f97316)
Glow: Orange shadow
Text: Black
Purpose: Show two-base polymorphisms
```

#### **TNPs/MNPs (Very Rare)**
```css
Background: Red (#dc2626)
Glow: Red shadow
Text: White
Purpose: Emphasize complex mutations
```

---

## 📏 Format Details

### Display Format
- **60 bases per line** (standard in bioinformatics)
- **10 lines per block** (600 bases visible at once)
- **Position markers** every 60 bases (1-indexed)
- **Match indicators** between sequences

### Why 60 Bases?
- Traditional alignment standard
- Easy to read and count
- Fits most screen widths
- Divisible by 3 (codons align nicely)

### Reading Frame Friendly
```
Position:  1                                                      60
Codon:     1       5        10        15        20
Twin 1:    ATG AAT TTA TCT GCT CTT CGC GTT GAA GAA GTA CAA AAT GTC...
           ||| ||| ||| ||| ||| ||| ||| ||| ||| ||| ||| ||| ||| |||
Twin 2:    ATG AAT TTA TCT GCT CTT CGC GTT GAA GAA GTA CAA AAT GTC...

Every 3 bases = 1 codon = 1 amino acid
```

---

## 💡 Use Cases

### 1. **Identifying Mutation Hotspots**
```
Use: Research disease susceptibility
Look for: Clusters of variants in density map
Action: Click cluster, examine variants
Result: Find critical mutation regions
```

### 2. **Conserved Region Analysis**
```
Use: Identify functionally important sequences
Look for: Long stretches of perfect matches
Interpretation: Conserved = probably critical
```

### 3. **Mutation Pattern Analysis**
```
Use: Understand mutation types
Look for: SNPs vs MNPs distribution
Pattern: Isolated SNPs = random mutations
         Clustered MNPs = recombination event
```

### 4. **Quality Control**
```
Use: Verify sequencing quality
Look for: Unusual patterns
Warning signs:
  - Too many variants (>5%)
  - Large MNPs (check for errors)
  - Regular periodic patterns (artifact?)
```

### 5. **Educational Purposes**
```
Use: Teaching genetics
Show: Real twin sequences side-by-side
Demonstrate: How similar/different twins are
Learn: What mutations look like in practice
```

---

## 🔬 Scientific Accuracy

### Alignment Algorithm
- **Simple position-by-position comparison**
- No gaps or insertions (assumes same length)
- Perfect for identical twins (same genome length)
- Fast and efficient for browser display

### Limitations
```
✅ Good for:
  - Twin comparisons
  - Same-length sequences
  - Point mutations (SNPs)
  - Small polymorphisms

❌ Not ideal for:
  - Large insertions/deletions
  - Different species
  - Sequences with many gaps
  - Structural variants
```

### For Advanced Analysis
If you need:
- Gap alignments → Use BLAST or EMBOSS
- Phylogenetic trees → Use MEGA or RAxML
- Multiple sequence alignment → Use Clustal Omega
- Structural variants → Use specialized tools

---

## 📈 Performance

### Optimization Features
```
✅ Displays up to 10,000 bases smoothly
✅ Chunked rendering (10 lines at a time)
✅ Lazy loading for long sequences
✅ Efficient event handling
✅ Hardware-accelerated CSS
```

### Recommended Sequence Lengths
- **< 1,000 bases**: Instant, perfect
- **1,000 - 5,000 bases**: Fast, very usable
- **5,000 - 10,000 bases**: Good, may scroll
- **> 10,000 bases**: Use detailed analysis instead

### Browser Requirements
- Modern browser (Chrome, Firefox, Safari, Edge)
- JavaScript enabled
- CSS3 support
- Decent screen resolution (1280x720+)

---

## 🎓 Reading an Alignment

### Step-by-Step Guide

#### 1. Orient Yourself
```
Check position numbers (left side)
Find start (position 1)
Note total length
```

#### 2. Scan for Variants
```
Look for colored highlights
Yellow/orange = SNPs (common)
Red = Complex mutations (rare)
```

#### 3. Check Match Indicators
```
|||||||  = Perfect match
||||×|||  = One difference
×××××××  = Many differences (check!)
```

#### 4. Use Density Map
```
Find high bars = variant hotspots
Click bar to navigate
Compare distribution
```

#### 5. Investigate Details
```
Click interesting positions
Read codon/frame info
Assess biological impact
```

---

## 🧬 Real Example Walkthrough

### Sample: BRCA1 Twin Comparison

#### What You'll See:
```
Position 1-600:
  - Mostly matches (||||||||)
  - Background light gray
  - Sequences look identical
  
Position 645-650:
  - Bright red highlight! 🔴
  - MNP detected
  - 6 bases different
  - Affects codons 215-216
  
Position 651-780:
  - Back to matches
  - Normal gray background
  - Sequences identical again
```

#### Interpretation:
```
✅ 99.2% identity overall
⚠️  One complex mutation at position 645-650
✅ Rest of gene perfectly conserved
🔬 Conclusion: One significant variant to investigate
```

---

## 💻 Keyboard Shortcuts

While in alignment view:

```
Arrow Keys:  Navigate sequence
Home/End:    Jump to start/end
Page Up/Down: Scroll by block
Esc:         Close tooltip
```

---

## 🎯 Pro Tips

### Tip 1: Use Density Map First
```
Before reading sequence:
1. Look at density map
2. Identify interesting regions
3. Click to jump there
4. Read in detail

Saves time on long sequences!
```

### Tip 2: Check Match Indicators
```
Quick quality check:
- Mostly ||| = Good alignment
- Many ××× = Check for errors
- Regular patterns = Possible artifact
```

### Tip 3: Screenshot for Records
```
Navigate to interesting region
Zoom browser to fit on screen
Take screenshot
Document findings
```

### Tip 4: Compare to Reference
```
If you have reference genome:
Load as Twin 1 (reference)
Load sample as Twin 2 (test)
Variants show differences from reference
```

### Tip 5: Look for Patterns
```
Clustered variants = Hotspot/recombination
Isolated SNPs = Random mutations
Regular spacing = Check for artifacts
No variants = Conserved/critical region
```

---

## 🆘 Troubleshooting

### Alignment Not Showing
**Check:**
- Did you load 2 sequences? (FASTA format)
- Are they similar lengths?
- Did analysis complete?

### Can't See Variants
**Try:**
- Scroll down to alignment section
- Check legend colors
- Look at density map for guidance
- Zoom browser if text too small

### Slow Performance
**Solutions:**
- Use shorter sequences (<5000 bp)
- Close other browser tabs
- Update browser to latest version
- Try desktop instead of mobile

### Tooltip Not Appearing
**Fix:**
- Make sure JavaScript enabled
- Try different browser
- Check browser console for errors
- Refresh page

---

## 📚 Further Reading

### Learn More About:
- **Sequence Alignment**: NCBI BLAST tutorial
- **Variant Calling**: SAMtools documentation
- **Twin Genetics**: TwinsUK research papers
- **BRCA1 Gene**: ClinVar database
- **Bioinformatics**: Rosalind.info exercises

### Tools for Advanced Analysis:
- **BLAST**: Sequence similarity search
- **Clustal Omega**: Multiple alignment
- **IGV**: Genome viewer
- **UCSC Genome Browser**: Comprehensive viewer

---

## ✅ Quick Checklist

**Before Using:**
- [ ] Load twin sequences (FASTA format)
- [ ] Run analysis
- [ ] Wait for results to appear

**While Exploring:**
- [ ] Check density map for overview
- [ ] Scan for colored highlights
- [ ] Hover over bases for details
- [ ] Click for position information
- [ ] Note interesting variants

**After Analysis:**
- [ ] Screenshot key regions
- [ ] Document variant positions
- [ ] Look up in databases (ClinVar, etc.)
- [ ] Assess biological significance
- [ ] Share findings!

---

## 🎉 Summary

The Full Sequence Alignment View gives you:

✅ **Complete picture** of both sequences side-by-side
✅ **Visual highlighting** of all variants in context  
✅ **Interactive exploration** with hover and click
✅ **Density map** showing variant distribution
✅ **Professional quality** alignment display
✅ **Fast navigation** for long sequences
✅ **Detailed tooltips** with position information
✅ **Color-coded** variant classification

**It's like having BLAST, IGV, and a sequence viewer all in one, right in your browser!** 🚀

---

**Try it now with the twin sample data and see your sequences come to life!** 🧬✨

