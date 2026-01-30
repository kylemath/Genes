# 🧬 DNA Analysis Playground

🚀 **[Live Demo](https://kylemath.github.io/Genes)** 🚀

An interactive web application for learning computational genomics and bioinformatics. No installation required - just open in your browser!

## 🚀 Quick Start

1. Open `index.html` in your web browser (just double-click the file)
2. Click one of the sample sequences or paste your own DNA sequence
3. Click "Analyze Sequence" to see the results!

## 📚 What You'll Learn

This tool demonstrates fundamental DNA analysis techniques:

### 1. **FASTA Format**
The standard format for representing DNA sequences, starting with a header line (>) followed by the sequence.

### 2. **Nucleotide Counting**
Basic composition analysis - counting A, T, G, and C bases in your sequence.

### 3. **GC Content**
The percentage of Guanine (G) and Cytosine (C) in the sequence. This is important because:
- GC pairs have 3 hydrogen bonds vs. 2 for AT pairs
- Affects DNA melting temperature
- Influences gene expression and chromosome structure
- Different organisms have characteristic GC contents

### 4. **Reverse Complement**
DNA is double-stranded with complementary base pairing (A↔T, G↔C). The reverse complement shows the sequence of the opposite strand read in the opposite direction.

### 5. **Open Reading Frames (ORFs)**
Regions that could potentially code for proteins:
- Start with ATG (start codon)
- End with TAA, TAG, or TGA (stop codons)
- Read in triplets (codons)

### 6. **Codon Analysis**
Three-nucleotide sequences that code for amino acids:
- 64 possible codons
- Code for 20 amino acids
- Include 3 stop signals

## 🧪 Sample Sequences Included

1. **E. coli LacZ Gene Fragment**: Part of the beta-galactosidase gene (famous in molecular biology!)
2. **Human Insulin Gene Fragment**: Part of the gene that codes for insulin
3. **Simple Test Sequence**: A short sequence for quick testing

## 🔬 Features

- **Interactive Visualizations**: Colorful charts and graphs
- **Real-time Analysis**: Instant results as you input sequences
- **Educational Explanations**: Learn what each analysis means
- **No Dependencies**: Pure HTML, CSS, and JavaScript
- **Mobile Friendly**: Works on phones and tablets

## 📖 How to Use Your Own Sequences

### Option 1: Paste DNA Sequence
Simply paste any DNA sequence (A, T, G, C letters) into the text box.

### Option 2: Load FASTA File
Click "Load FASTA File" and select a `.fasta` or `.fa` file from your computer.

### Where to Get More Sequences

1. **NCBI GenBank**: https://www.ncbi.nlm.nih.gov/genbank/
   - The world's largest public DNA database
   - Search for any gene or organism

2. **Ensembl**: https://www.ensembl.org/
   - Genome browser for vertebrates and other species

3. **UCSC Genome Browser**: https://genome.ucsc.edu/
   - Extensive collection of genome sequences

## 🎓 Learning Path

### Beginner
1. Start with the "Simple Test Sequence"
2. Understand nucleotide counting and GC content
3. Explore the reverse complement

### Intermediate
1. Try the E. coli or insulin samples
2. Understand Open Reading Frames
3. Learn about codon usage

### Advanced
1. Download real genome sequences from NCBI
2. Compare GC content across different organisms
3. Find potential genes using ORF analysis

## 🛠️ Technical Details

### Technologies Used
- **HTML5**: Structure and semantic markup
- **CSS3**: Modern styling with gradients and animations
- **Vanilla JavaScript**: No frameworks or libraries needed

### Browser Compatibility
Works in all modern browsers:
- Chrome/Edge (recommended)
- Firefox
- Safari
- Opera

### File Structure
```
Genes/
├── index.html      # Main HTML structure
├── styles.css      # All styling and visual design
├── script.js       # Analysis algorithms and logic
└── README.md       # This file
```

## 🔍 Understanding the Results

### Nucleotide Composition Chart
Shows the count of each base in your sequence. Useful for:
- Detecting sequencing bias
- Understanding sequence complexity

### GC Content Pie Chart
Visualizes the ratio of GC to AT bases. Important for:
- Identifying gene-rich regions (often higher GC)
- PCR primer design
- Predicting DNA stability

### Open Reading Frames (ORFs)
Lists potential protein-coding regions:
- **Position**: Where in the sequence it starts and ends
- **Length**: Size in base pairs and codons
- **Reading Frame**: Which of 3 possible frames (+1, +2, +3)
- **Protein**: Translated amino acid sequence

### Codon Frequency
Shows which three-letter codes are most common:
- Different organisms prefer different codons
- Useful for optimizing genes for expression
- Can indicate evolutionary relationships

## 🧬 Fun Facts

- The human genome contains about 3 billion base pairs
- If you uncoiled all DNA in your body, it would stretch to the Sun and back multiple times
- E. coli has about 4.6 million base pairs
- The smallest known genome is about 160,000 base pairs (certain bacteria)

## 📚 Further Learning

### Books
- "Molecular Biology of the Cell" by Alberts et al.
- "Bioinformatics for Dummies" by Jean-Michel Claverie

### Online Courses
- Coursera: "Biology Meets Programming: Bioinformatics for Beginners"
- edX: "Introduction to Genomic Technologies"

### Tools & Databases
- **BLAST**: Compare sequences (https://blast.ncbi.nlm.nih.gov/)
- **UniProt**: Protein sequence database
- **PDB**: Protein structure database

## 🤝 Contributing

Feel free to extend this project! Ideas:
- Add more sample sequences
- Implement additional analyses (e.g., restriction sites, melting temperature)
- Create protein structure predictions
- Add sequence alignment tools

## 📄 License

This is an educational project. Feel free to use, modify, and share!

## 🙋 Need Help?

### Common Issues

**Q: Why don't I see any ORFs?**
A: Your sequence might be too short or not contain complete ORFs (need ATG start and stop codon). Try one of the sample sequences.

**Q: Can I analyze RNA sequences?**
A: This tool is designed for DNA. For RNA, replace T with U in your sequence first.

**Q: What's the maximum sequence length?**
A: The tool can handle sequences up to several megabases, but visualization is limited to the first 500 nucleotides for performance.

---

**Happy analyzing! 🧬🔬**

