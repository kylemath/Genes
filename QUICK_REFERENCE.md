# 🧬 DNA Analysis Quick Reference Card

## 📁 Your Project Files

| File | Purpose |
|------|---------|
| `index.html` | Main DNA analyzer web app |
| `data-access.html` | Portal to access real twin data |
| `styles.css` | Visual styling |
| `script.js` | Analysis algorithms |
| `sample_twin_data.fasta` | Example twin sequences |
| `DATA_SOURCES_GUIDE.md` | Comprehensive data access guide |
| `TWIN_ANALYSIS_TUTORIAL.md` | Twin study tutorial |
| `README.md` | Main documentation |

---

## 🚀 Quick Start Commands

### Open the App
```bash
# Just double-click index.html in Finder
# OR navigate in terminal:
cd ~/Coding/Genes
open index.html
```

### View Data Portal
```bash
open data-access.html
```

---

## 🌐 Essential Websites (Bookmarks)

### Public Databases
- **GEO**: https://www.ncbi.nlm.nih.gov/geo/ (Gene expression data)
- **SRA**: https://www.ncbi.nlm.nih.gov/sra (Raw sequencing)
- **NCBI Gene**: https://www.ncbi.nlm.nih.gov/gene/ (Gene info)
- **1000 Genomes**: https://www.internationalgenome.org/ (Population data)
- **Ensembl**: https://www.ensembl.org/ (Genome browser)

### Analysis Tools
- **BLAST**: https://blast.ncbi.nlm.nih.gov/ (Compare sequences)
- **EMBOSS**: https://www.ebi.ac.uk/Tools/emboss/ (Sequence analysis)
- **Galaxy**: https://usegalaxy.org/ (Workflow platform)

### Twin Data Sources
- **TwinsUK**: https://twinsuk.ac.uk/
- **Australian Twins**: https://www.twins.org.au/
- **dbGaP**: https://www.ncbi.nlm.nih.gov/gap/ (Controlled access)

---

## 📊 Key Twin Datasets

| Dataset ID | Description | Difficulty | Link |
|------------|-------------|------------|------|
| GSE100227 | Australian twins (132 pairs) | Easy | [Access](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100227) |
| GSE65638 | Chinese twins (8 subjects) | Easy | [Access](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65638) |
| GSE61862 | Schizophrenia twins | Medium | [Access](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61862) |

---

## 🔬 DNA Analysis Cheat Sheet

### Basic Terms
- **Nucleotide**: A, T, G, C (DNA building blocks)
- **Codon**: 3 nucleotides = 1 amino acid
- **ORF**: Open Reading Frame (potential gene)
- **GC Content**: % of G and C bases (important for stability)

### File Formats

#### FASTA (Sequence)
```
>Sequence_Name
ATGCGTACGTAGCTAGCT
```

#### FASTQ (Sequencing)
```
@Read_ID
ATGCGTACGT
+
!''*((((***
```

#### VCF (Variants)
```
#CHROM  POS  ID      REF  ALT
chr1    1234 rs123   A    G
```

### Genetic Code Quick Reference

| Codon | Amino Acid | Codon | Amino Acid |
|-------|------------|-------|------------|
| ATG   | Met (START)| TAA   | STOP       |
| TGG   | Trp        | TAG   | STOP       |
| AAA   | Lys        | TGA   | STOP       |
| GCT   | Ala        | TTT   | Phe        |

---

## 🛠️ Command Line Tools (Optional)

### Installation (macOS)
```bash
# Install Homebrew first
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install tools
brew install sra-tools  # Download sequences
brew install emboss     # Analysis tools
brew install samtools   # BAM/SAM files
brew install seqtk      # FASTQ/FASTA conversion
```

### Basic Commands

```bash
# Count sequences in FASTA
grep -c "^>" file.fasta

# Convert FASTQ to FASTA
seqtk seq -a input.fastq > output.fasta

# Download from SRA
fastq-dump --split-files SRR1234567

# Align two sequences
needle seq1.fasta seq2.fasta -outfile alignment.txt

# Extract specific region
samtools faidx genome.fasta chr1:1000-2000
```

---

## 📖 Analysis Workflow

### For Beginners (Web App Only)
1. Open `index.html`
2. Click a sample (try "Twin BRCA1 Comparison")
3. Click "Analyze"
4. Explore results
5. Read explanations

### For Intermediate (Real Data)
1. Visit `data-access.html`
2. Click on a dataset (GSE100227 recommended)
3. Download supplementary files
4. Open in Excel/text editor
5. Find sequence data
6. Paste into web app

### For Advanced (Command Line)
1. Find dataset in SRA
2. Download with `fastq-dump`
3. Convert to FASTA with `seqtk`
4. Analyze with tools
5. Visualize in web app

---

## 🎓 Learning Path

### Week 1: Basics
- [ ] Complete web app tutorial
- [ ] Understand nucleotides, codons, ORFs
- [ ] Analyze all sample sequences
- [ ] Read README.md

### Week 2: Real Data
- [ ] Download first GEO dataset
- [ ] Explore in spreadsheet
- [ ] Read DATA_SOURCES_GUIDE.md
- [ ] Find gene sequences

### Week 3: Twin Studies
- [ ] Read TWIN_ANALYSIS_TUTORIAL.md
- [ ] Compare twin sequences
- [ ] Calculate similarities
- [ ] Read 1 research paper

### Week 4: Advanced
- [ ] Install command-line tools
- [ ] Download SRA data
- [ ] Perform alignments
- [ ] Start your own project!

---

## 🆘 Troubleshooting

### Problem: Web app doesn't work
**Solution**: Make sure all files are in same folder, open in modern browser (Chrome/Firefox)

### Problem: Can't download dataset
**Solution**: Some require registration, try different dataset or use browser's "Save As"

### Problem: File too large
**Solution**: Start with smaller regions, use command-line tools to extract portions

### Problem: Don't understand results
**Solution**: Check explanations in web app, read README.md, ask on Biostars

---

## 💡 Analysis Tips

### Comparing Sequences
1. Always align before comparing
2. Note sequence lengths
3. Check quality scores (FASTQ)
4. Consider multiple alignments

### Finding Mutations
1. Use reference sequence
2. Look for SNVs (single changes)
3. Check if synonymous or not
4. Verify it's not sequencing error

### Twin Analysis
1. Start with identical twins
2. Compare same tissue/cell type
3. Control for age
4. Consider environment

---

## 📱 Mobile Access

The web app works on mobile! Just:
1. Transfer files to phone
2. Open in mobile browser
3. Smaller datasets recommended
4. Use landscape orientation

---

## 🔗 Useful Search Queries

### PubMed
```
"twin study" AND genetics AND [disease]
epigenetics AND monozygotic twins
heritability AND twin cohort
```

### Google Scholar
```
twin concordance [trait]
identical twins genetic differences
twin registry genomics
```

### NCBI GEO
```
twins methylation
monozygotic dizygotic
twin discordant
```

---

## 📊 Statistical Formulas

### GC Content
```
GC% = ((G + C) / (A + T + G + C)) × 100
```

### Sequence Identity
```
Identity% = (Matches / Alignment_Length) × 100
```

### Heritability (Falconer's)
```
h² = 2 × (Concordance_MZ - Concordance_DZ)
```

### Concordance Rate
```
Concordance = Both_Affected / At_Least_One_Affected
```

---

## 🎯 Project Ideas

### Beginner Projects
1. Compare GC content across species
2. Find longest ORF in a gene
3. Translate DNA to protein
4. Count codon usage

### Intermediate Projects
1. Analyze disease-associated genes
2. Compare twin methylation data
3. Find conserved regions
4. Study mutation effects

### Advanced Projects
1. Whole genome analysis
2. Population genomics
3. Phylogenetic analysis
4. Machine learning on sequences

---

## 📚 Recommended Reading Order

1. **README.md** - Start here!
2. **Web app** - Hands-on learning
3. **DATA_SOURCES_GUIDE.md** - Find data
4. **TWIN_ANALYSIS_TUTORIAL.md** - Deep dive
5. **Research papers** - Current science

---

## 🌟 Best Practices

### Data Management
- Keep original files
- Document everything
- Use version control
- Organize by project

### Analysis
- Start simple
- Validate results
- Check for errors
- Replicate findings

### Learning
- Read documentation
- Try examples first
- Ask questions
- Join communities

---

## 📞 Getting Help

### Forums
- **Biostars**: https://www.biostars.org/ (Best for bioinformatics)
- **Reddit**: r/bioinformatics, r/genetics
- **Stack Overflow**: For programming questions

### Documentation
- NCBI Handbook: https://www.ncbi.nlm.nih.gov/books/
- Ensembl Help: https://www.ensembl.org/Help/
- Galaxy Training: https://training.galaxyproject.org/

---

## ✅ Quick Checklist

**Getting Started:**
- [ ] Open index.html
- [ ] Try all samples
- [ ] Read explanations
- [ ] Check data-access.html

**First Analysis:**
- [ ] Choose dataset
- [ ] Download data
- [ ] Extract sequences
- [ ] Analyze & interpret

**Going Further:**
- [ ] Install tools
- [ ] Join forum
- [ ] Read papers
- [ ] Start project

---

## 🎓 Key Concepts to Master

1. **DNA Structure**: Double helix, base pairing, antiparallel
2. **Gene Expression**: DNA → RNA → Protein
3. **Mutations**: Types and effects
4. **Twin Studies**: MZ vs DZ, concordance, heritability
5. **Epigenetics**: Beyond the sequence
6. **Bioinformatics**: Algorithms and tools

---

**Print this page for quick reference! 📋**

**Questions? Check the guides or visit Biostars.org! 🚀**

