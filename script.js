// DNA Analysis Playground - Main JavaScript File

// Sample DNA sequences for educational purposes
const samples = {
    ecoli: {
        name: "E. coli LacZ Gene Fragment",
        description: "Fragment from the beta-galactosidase gene",
        sequence: `>E.coli_LacZ_fragment
ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCT
GGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGC
GAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGC
TTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCT
GAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATC
TACACCAACGTGACCTATCCCATTACGGTCAATCCGCCGTTTGTTCCCACGGAGAATCCG
ACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCTGGCTACAGGAAGGCCAGACG
CGAATTATTTTTGATGGCGTTAACTCGGCGTTTCATCTGTGGTGCAACGGGCGCTGGGTC
GGTTACGGCCAGGACAGTCGTTTGCCGTCTGAATTTGACCTGAGCGCATTTTTACGCGCC
GGAGAAAACCGCCTCGCGGTGATGGTGCTGCGTTGGAGTGACGGCAGTTATCTGGAAGAT
CAGGATATGTGGCGGATGAGCGGCATTTTCCGTGACGTCTCGTTGCTGCATAAACCGACT
ACACAAATCAGCGATTTCCATGTTGCCACTCGCTTTAATGATGATTTCAGCCGCGCTGTA
CTGGAGGCTGAAGTTCAGATGTGCGGCGAGTTGCGTGACTACCTACGGGTAACAGTTTCT
TTATGGCAGGGTGAAACGCAGGTCGCCAGCGGCACCGCGCCTTTCGGCGGTGAAATTATC
GATGAGCGTGGTGGTTATGCCGATCGCGTCACACTACGTCTGAACGTCGAAAACCCGAAA
CTGTGGAGCGCCGAAATCCCGAATCTCTATCGTGCGGTGGTTGAACTGCACACCGCCGAC
GGCACGCTGATTGAAGCAGAAGCCTGCGATGTCGGTTTCCGCGAGGTGCGGATTGAAAAT
GGTCTGCTGCTGCTGAACGGCAAGCCGTTGCTGATTCGAGGCGTTAACCGTCACGAGCAT
CATCCTCTGCATGGTCAGGTCATGGATGAGCAGACGATGGTGCAGGATATCCTGCTGATG
AAGCAGAACAACTTTAACGCCGTGCGCTGTTCGCATTATCCGAACCATCCGCTGTGGTAC
ACGCTGTGCGACCGCTACGGCCTGTATGTGGTGGATGAAGCCAATATTGAAACCCACGGC
ATGGTGCCAATGAATCGTCTGACCGATGATCCGCGCTGGCTACCGGCGATGAGCGAACGC
GTAACGCGAATGGTGCAGCGCGATCGTAATCACCCGAGTGTGATCATCTGGTCGCTGGGG
AATGAATCAGGCCACGGCGCTAATCACGACGCGCTGTATCGCTGGATCAAATCTGTCGAT
CCTTCCCGCCCGGTGCAGTATGAAGGCGGCGGAGCCGACACCACGGCCACCGATATTATT
TGCCCGATGTACGCGCGCGTGGATGAAGACCAGCCCTTCCCGGCTGTGCCGAAATGGTCC
ATCAAAAAATGGCTTTCGCTACCTGGAGAGACGCGCCCGCTGATCCTTTGCGAATACGCC
CACGCGATGGGTAACAGTCTTGGCGGTTTCGCTAAATACTGGCAGGCGTTTCGTCAGTAT
CCCCGTTTACAGGGCGGCTTCGTCTGGGACTGGGTGGATCAGTCGCTGATTAAATATGAT
GAAAACGGCAACCCGTGGTCGGCTTACGGCGGTGATTTTGGCGATACGCCGAACGATCGC
CAGTTCTGTATGAACGGTCTGGTCTTTGCCGACCGCACGCCGCATCCAGCGCTGACGGAA
GCAAAACACCAGCAGCAGTTTTTCCAGTTCCGTTTATCCGGGCAAACCATCGAAGTGACC
AGCGAATACCTGTTCCGTCATAGCGATAACGAGCTCCTGCACTGGATGGTGGCGCTGGAT
GGTAAGCCGCTGGCAAGCGGTGAAGTGCCTCTGGATGTCGCTCCACAAGGTAAACAGTTG
ATTGAACTGCCTGAACTACCGCAGCCGGAGAGCGCCGGGCAACTCTGGCTCACAGTACGC
GTAGTGCAACCGAACGCGACCGCATGGTCAGAAGCCGGGCACATCAGCGCCTGGCAGCAG
TGGCGTCTGGCGGAAAACCTCAGTGTGACGCTCCCCGCCGCGTCCCACGCCATCCCGCAT
CTGACCACCAGCGAAATGGATTTTTGCATCGAGCTGGGTAATAAGCGTTGGCAATTTAAC
CGCCAGTCAGGCTTTCTTTCACAGATGTGGATTGGCGATAAAAAACAACTGCTGACGCCG
CTGCGCGATCAGTTCACCCGTGCACCGCTGGATAACGACATTGGCGTAAGTGAAGCGACC
CGCATTGACCCTAACGCCTGGGTCGAACGCTGGAAGGCGGCGGGCCATTACCAGGCCGAA
GCAGCGTTGTTGCAGTGCACGGCAGATACACTTGCTGATGCGGTGCTGATTACGACCGCT
CACGCGTGGCAGCATCAGGGGAAAACCTTATTTATCAGCCGGAAAACCTACCGGATTGAT
GGTAGTGGTCAAATGGCGATTACCGTTGATGTTGAAGTGGCGAGCGATACACCGCATCCG
GCGCGGATTGGCCTGAACTGCCAGCTGGCGCAGGTAGCAGAGCGGGTAAACTGGCTCGGA
TTAGGGCCGCAAGAAAACTATCCCGACCGCCTTACTGCCGCCTGTTTTGACCGCTGGGAT
CTGCCATTGTCAGACATGTATACCCCGTACGTCTTCCCGAGCGAAAACGGTCTGCGCTGC
GGGACGCGCGAATTGAATTATGGCCCACACCAGTGGCGCGGCGACTTCCAGTTCAACATC
AGCCGCTACAGTCAACAGCAACTGATGGAAACCAGCCATCGCCATCTGCTGCACGCGGAA
GAAGGCACATGGCTGAATATCGACGGTTTCCATATGGGGATTGGTGGCGACGACTCCTGG
AGCCCGTCAGTATCGGCGGAATTCCAGCTGAGCGCCGGTCGCTACCATTACCAGTTGGTC
TGGTGTCAAAAATAATAATAACCGGGCAGGCCATGTCTGCCCGTATTTCGCGTAAGGAAA
TCCATTATGTACTATTTAAAAAACACAAACTTTTGGATGTTCGGTTTATTCTTTTTCTTTT
ACTTTTTTATCATGGGAGCCTACTTCCCGTTTTTCCCGATTTGGCTACATGACATCAACC
ATATCAGCAAAAGTGATACGGGTATTTTTGCCGCATTCTATCAGATTTATCAGGCAATAA
ATATTTATTATAAGTTGTCTGCTGCGCTGTGCTGATTGCGGCGATGTCGATTCTCACCTG
TACTTGCGAGCCGCTGGATAAGGTCAACTCGTTTCACCACTCCGAGAACGCCGGCCTCGA
TGCGTAACACGGTCGGTTATGGTTAAAATTAAGCACGACGGGAGAATATATGAGTGGCAA
TCTATGAGTGGTAGTTAATGATGTTTGAAACGAGTCGAATCTTATCGGCTACATCTCTCG
AACGTTTAATGGCAAATCTAAAACAATTGTCTTCCTTATAAAACCTGAAAATCGACATAG
ATGGCAGATGCTAACCCGAGAGGGACTGTTGCTGATTTATGTGGCGGTTAGCGATAGTAA
TACGGATATATGTTAACTCGACTAACATATATGTGGAAGAGACTGACTTTTTAATACGAG
TCATCGGCGACAAGCCAGGAGAAACCTGTTGCGGATACGGGAAACGCATGAGTAAATTCT
TGAGTGAGAGCGGTAGGGCAACTTAGCGACAGTATTCAGGGGTTATTAGGAGGGATTATG
GCCTCTAAACGGGTCTTGAGGGGTTTTTTGCTGAAAGGAGGAACTATATCCGGAT`
    },
    insulin: {
        name: "Human Insulin Gene Fragment",
        description: "Part of the human insulin precursor gene",
        sequence: `>Human_Insulin_Precursor_Fragment
ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGAC
CCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTAC
CTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGAC
CTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTG
GCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGC
TCCCTCTACCAGCTGGAGAACTACTGCAACTAGACGCAGCCCGCAGGCAGCCCCACACCC
GCCGCCTCCTGCACCGAGAGAGATGGAATAAAGCCCTTGAACCAGC`
    },
    small: {
        name: "Simple Test Sequence",
        description: "A short sequence for quick testing",
        sequence: `>Test_Sequence
ATGGCTAGCGGTAAATCGTACGCTAGCTAAGGATCCGAATTCTGA`
    },
    twins: {
        name: "Twin Comparison Study",
        description: "Real BRCA1 gene fragments showing various mutation types",
        sequence: `>Twin_A_BRCA1_Fragment | Female | Age_32 | Healthy
ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAA
ATCTTAGAGTGTCCCATCTGGTAAGTCAGCACAGAGGGCCAGCTCATTCCTGCCCTCTGC
ACCGAGAGAGATGGAATAAAGCCCTTGAACCAGCCTTTGCAGAATGGAAACGTGCAGTCA
GCTGAGTGGCTGACAGAGCAGGGCTCCACACTGAGCTGTACCACATAGCCAGTGATGTAA
GTGAACAGGGAGAAACTGACCTGGCAGAACCTGATGCTGTGGTGAATACAGATCCTGAAA
CTGGGCTCAGTGGATGGGAAGGGCAGTGGGAACAGATGATGGGTCTACAGTTTCTGGGGA
AGTGAAGACAGAAATGCTAACCCCAAATGGGGCAAAATAAAGCAATTGCATGAACCAGAG
AGGGCCTCCTTCCAGTTTCCCTGGCTCTGGGCCCTCTGGTTTTCCATCTCTCTCCATCCC
CGGGAACTAGAAGTTGGCAGTTCAGACACTCAGCTCAGGGGAAGCAGGTGGAAAGGTCCC
GGGGGCTCCAGAGTGAGGCCCCTAGAATTTTGCTGGCAGAGTCTCAATTGGGAAGGCAGA
AGAGAATCTCTGATCTCTCTGTTAGAAATAGGTGTTATAAAAGCAGCATTTGGAAGATGC
CTATGCAAAAGAATCTCCCCAAGATCCTCTGTGGCTACAGATTGAGAAGAAGAGGCTGTG
CTCGAGTTTTGAAAGCAAAAATGGTCAAACAAGGTGACACTTAGAAGGAAAAACGGCTTG

>Twin_B_BRCA1_Fragment | Female | Age_32 | Breast_Cancer_Diagnosed
ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAA
ATCTTAGAGTGTCCCATCTGGTAAGTCAGCACAGAGGGCCAGCTCATTCCTGCCCTCTGC
ACCGAGAGAGATGGAATAAAGCCCTTGAACCAGCCTTTGCAGAATGGAAACGTGCAGTCA
GCTGAGTGGCTGACAGAGCAGGGCTCCACACTGAGCTGTACCACATAGCCAGTGATGTAA
GTGAACAGGGAGAAACTGACCTGGCAGAACCTGATGCTGTGGTGAATACAGATCCTGAAA
CTGGGCTCAGTGGATGGGAAGGGCAGTGGGAACAGATGATGGGTCTACAGTTTCTGGGGA
AGTGAAGACAGAAATGCTAACCCCAAATGGGGCAAAATAAAGCAATTGCATGAACCAGAG
AGGGCCTCCTTCCAGTTTCCCTGGCTCTGGGCCCTCTGGTTTTCCATCTCTCTCCATCCC
CGGGAACTAGAAGTTGGCAGTTCAGACACTCAGCTCAGGGGAAGCAGGTGGAAAGGTCCC
GGGGGCTCCAGAGTGAGGCCCCTAGAATTTTGCTGGCAGAGTCTCAATTGGGAAGGCAGA
AGAGAATCTCTGATCTCTCTGTTAGAAATAGGTGTTATAAAAGCAGCATTTGGGCTTAAC
CTATGCAAAAGAATCTCCCCAAGATCCTCTGTGGCTACAGATTGAGAAGAAGAGGCTGTG
CTCGAGTTTTGAAAGCAAAAATGGTCAAACAAGGTGACACTTAGAAGGAAAAACGGCTTG`
    }
};

// Known genes on chromosome 21 (major genes with clinical significance)
const chr21Genes = [
    { name: 'APP', start: 25880550, end: 26171128, description: 'Amyloid Beta Precursor Protein (Alzheimer\'s)' },
    { name: 'SOD1', start: 31659622, end: 31668931, description: 'Superoxide Dismutase (ALS)' },
    { name: 'DYRK1A', start: 37365666, end: 37530949, description: 'Dual-Specificity Kinase (Down Syndrome)' },
    { name: 'CSTB', start: 43776148, end: 43779192, description: 'Cystatin B (Epilepsy)' },
    { name: 'RUNX1', start: 34787801, end: 35062446, description: 'Runt-Related Transcription Factor (Leukemia)' },
    { name: 'CBS', start: 43354370, end: 43380062, description: 'Cystathionine Beta-Synthase (Homocystinuria)' },
    { name: 'CXADR', start: 17496471, end: 17532607, description: 'Coxsackievirus and Adenovirus Receptor' },
    { name: 'ITSN1', start: 33321483, end: 33530401, description: 'Intersectin 1 (Down Syndrome)' },
    { name: 'COL6A1', start: 45971413, end: 46015284, description: 'Collagen Type VI Alpha 1' },
    { name: 'COL6A2', start: 46017395, end: 46041235, description: 'Collagen Type VI Alpha 2' },
    { name: 'S100B', start: 47405643, end: 47412782, description: 'Calcium Binding Protein (Down Syndrome marker)' },
    { name: 'GART', start: 33243058, end: 33271569, description: 'Phosphoribosylglycinamide Formyltransferase' },
    { name: 'SON', start: 33032978, end: 33072358, description: 'SON DNA Binding Protein' },
    { name: 'DOPEY2', start: 37087778, end: 37181618, description: 'Dopey Family Member 2' },
    { name: 'KCNJ6', start: 38776577, end: 39180820, description: 'Potassium Channel (Neuronal)' },
    { name: 'TMPRSS2', start: 41464308, end: 41507950, description: 'Transmembrane Serine Protease' },
    { name: 'ERG', start: 38380525, end: 38670044, description: 'ETS Transcription Factor' },
    { name: 'ETS2', start: 38780149, end: 38860043, description: 'ETS Transcription Factor 2' },
    { name: 'HMGN1', start: 38435620, end: 38441489, description: 'High Mobility Group Nucleosomal Binding Domain 1' },
    { name: 'AIRE', start: 44511530, end: 44524124, description: 'Autoimmune Regulator' }
];

// Find which gene a position falls in
function findGeneAtPosition(position, chromosome = '21') {
    if (chromosome !== '21' && chromosome !== 'chr21') {
        return null;
    }
    
    for (const gene of chr21Genes) {
        if (position >= gene.start && position <= gene.end) {
            return gene;
        }
    }
    return null;
}

// Genetic Code Table (Standard)
const geneticCode = {
    'TTT': 'Phe', 'TTC': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
    'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser',
    'TAT': 'Tyr', 'TAC': 'Tyr', 'TAA': 'STOP', 'TAG': 'STOP',
    'TGT': 'Cys', 'TGC': 'Cys', 'TGA': 'STOP', 'TGG': 'Trp',
    'CTT': 'Leu', 'CTC': 'Leu', 'CTA': 'Leu', 'CTG': 'Leu',
    'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAT': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met',
    'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAT': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val',
    'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
};

// File input handler with large file support
document.getElementById('fileInput').addEventListener('change', function(event) {
    const file = event.target.files[0];
    if (!file) return;

    const fileSizeMB = file.size / (1024 * 1024);
    
    // For files larger than 10MB, parse directly without showing in textarea
    if (fileSizeMB > 10) {
        const reader = new FileReader();
        reader.onload = function(e) {
            const text = e.target.result;
            
            // Don't put large text in textarea - parse directly
            document.getElementById('sequenceInput').value = `[Large file loaded: ${file.name} (${fileSizeMB.toFixed(1)} MB)]\nParsing... Please wait.`;
            
            // Parse and analyze directly
            setTimeout(() => {
                analyzeLargeFile(text);
            }, 100);
        };
        
        // Show loading message
        document.getElementById('sequenceInput').value = `Loading large file: ${file.name} (${fileSizeMB.toFixed(1)} MB)...\nThis may take a moment.`;
        
        reader.readAsText(file);
    } else {
        // Small file - normal behavior
        const reader = new FileReader();
        reader.onload = function(e) {
            const text = e.target.result;
            document.getElementById('sequenceInput').value = text;
        };
        reader.readAsText(file);
    }
});

// Analyze large file directly without textarea
function analyzeLargeFile(fileContent) {
    try {
        // Parse FASTA
        const parsed = parseFasta(fileContent);
        
        // Update textarea with summary instead of full content
        if (Array.isArray(parsed) && parsed.length >= 2) {
            const seq1Len = parsed[0].sequence.length;
            const seq2Len = parsed[1].sequence.length;
            document.getElementById('sequenceInput').value = 
                `✅ Large file parsed successfully!\n\n` +
                `Sequence 1: ${parsed[0].header}\nLength: ${seq1Len.toLocaleString()} bases\n\n` +
                `Sequence 2: ${parsed[1].header}\nLength: ${seq2Len.toLocaleString()} bases\n\n` +
                `Analyzing...`;
        }
        
        // Run analysis
        if (Array.isArray(parsed) && parsed.length >= 2) {
            analyzeTwinComparison(parsed);
        } else {
            const { header, sequence } = Array.isArray(parsed) ? parsed[0] : parsed;
            if (!sequence || sequence.length === 0) {
                alert('No valid DNA sequence found.');
                return;
            }
            // Single sequence analysis (existing code)
            analyzeSingleSequenceFromLargeFile(header, sequence);
        }
        
        // Update textarea with completion message
        setTimeout(() => {
            document.getElementById('sequenceInput').value = 
                `✅ Analysis complete!\n\n` +
                `Check results below ⬇️\n\n` +
                `(Full sequence not shown to prevent browser lag)`;
        }, 500);
        
    } catch (error) {
        console.error('Error analyzing large file:', error);
        document.getElementById('sequenceInput').value = 
            `❌ Error parsing file: ${error.message}\n\nPlease ensure file is in valid FASTA format.`;
    }
}

// Single sequence analysis for large files
function analyzeSingleSequenceFromLargeFile(header, sequence) {
    // Perform analyses
    const counts = countNucleotides(sequence);
    const gcContent = calculateGC(sequence);
    const revComp = reverseComplement(sequence);
    const orfs = findORFs(sequence);
    const codonCounts = countCodons(sequence);
    
    // Display basic statistics
    const statsHtml = `
        <div class="stat-item">
            <span class="stat-value">${sequence.length.toLocaleString()}</span>
            <span class="stat-label">Base Pairs</span>
        </div>
        <div class="stat-item">
            <span class="stat-value" style="color: var(--nucleotide-a)">${counts.A.toLocaleString()}</span>
            <span class="stat-label">Adenine (A)</span>
        </div>
        <div class="stat-item">
            <span class="stat-value" style="color: var(--nucleotide-t)">${counts.T.toLocaleString()}</span>
            <span class="stat-label">Thymine (T)</span>
        </div>
        <div class="stat-item">
            <span class="stat-value" style="color: var(--nucleotide-g)">${counts.G.toLocaleString()}</span>
            <span class="stat-label">Guanine (G)</span>
        </div>
        <div class="stat-item">
            <span class="stat-value" style="color: var(--nucleotide-c)">${counts.C.toLocaleString()}</span>
            <span class="stat-label">Cytosine (C)</span>
        </div>
        <div class="stat-item">
            <span class="stat-value">${gcContent.toFixed(2)}%</span>
            <span class="stat-label">GC Content</span>
        </div>
    `;
    document.getElementById('basicStats').innerHTML = statsHtml;
    
    // Render charts
    renderBarChart(counts, 'nucleotideChart', 'nucleotide');
    renderGCChart(gcContent);
    
    // Display reverse complement (limited for large sequences)
    const showLength = Math.min(200, sequence.length);
    document.getElementById('reverseComplement').innerHTML = `
        <h3>Original Sequence (5' → 3'):</h3>
        <div style="margin: 10px 0; word-break: break-all;">${sequence.substring(0, showLength)}${sequence.length > showLength ? '...' : ''}</div>
        <h3>Reverse Complement (3' → 5'):</h3>
        <div style="margin: 10px 0; word-break: break-all;">${revComp.substring(0, showLength)}${revComp.length > showLength ? '...' : ''}</div>
        ${sequence.length > showLength ? `<p style="color: var(--text-secondary); text-align: center; margin-top: 10px;">Showing first ${showLength} of ${sequence.length.toLocaleString()} bases</p>` : ''}
    `;
    
    // Display ORFs (limited to first 10)
    if (orfs.length > 0) {
        let orfHtml = '<div class="orf-list">';
        orfs.slice(0, 10).forEach((orf, index) => {
            orfHtml += `
                <div class="orf-item">
                    <div class="orf-header">
                        <div class="orf-title">ORF ${index + 1}</div>
                    </div>
                    <div class="orf-details">
                        <div class="orf-detail"><strong>Position:</strong> ${orf.start.toLocaleString()} - ${orf.end.toLocaleString()}</div>
                        <div class="orf-detail"><strong>Length:</strong> ${orf.length.toLocaleString()} bp (${Math.floor(orf.length/3)} codons)</div>
                        <div class="orf-detail"><strong>Reading Frame:</strong> +${orf.frame}</div>
                    </div>
                    <div class="orf-sequence">
                        <strong>DNA:</strong> ${orf.sequence.substring(0, 90)}${orf.sequence.length > 90 ? '...' : ''}<br>
                        <strong>Protein:</strong> ${orf.protein}
                    </div>
                </div>
            `;
        });
        orfHtml += '</div>';
        
        if (orfs.length > 10) {
            orfHtml += `<p class="text-center mt-20" style="color: var(--text-secondary);">
                Showing 10 of ${orfs.length.toLocaleString()} ORFs found
            </p>`;
        }
        
        document.getElementById('orfResults').innerHTML = orfHtml;
    } else {
        document.getElementById('orfResults').innerHTML = `
            <p style="color: var(--text-secondary); text-align: center; padding: 20px;">
                No significant ORFs found (minimum length: 30 bp)
            </p>
        `;
    }
    
    // Display codon analysis (top 20)
    const sortedCodons = Object.entries(codonCounts)
        .sort((a, b) => b[1] - a[1])
        .slice(0, 20);
    
    if (sortedCodons.length > 0) {
        let codonHtml = '<div class="codon-grid">';
        sortedCodons.forEach(([codon, count]) => {
            const aa = geneticCode[codon] || 'Unknown';
            codonHtml += `
                <div class="codon-item">
                    <div class="codon-name">${codon}</div>
                    <div class="codon-count">${aa}</div>
                    <div class="codon-count">Count: ${count.toLocaleString()}</div>
                </div>
            `;
        });
        codonHtml += '</div>';
        
        if (Object.keys(codonCounts).length > 20) {
            codonHtml += `<p class="text-center mt-20" style="color: var(--text-secondary);">
                Showing top 20 of ${Object.keys(codonCounts).length} codons found
            </p>`;
        }
        
        document.getElementById('codonAnalysis').innerHTML = codonHtml;
    } else {
        document.getElementById('codonAnalysis').innerHTML = `
            <p style="color: var(--text-secondary); text-align: center; padding: 20px;">
                Sequence too short for codon analysis
            </p>
        `;
    }
    
    // Render sequence visualization (optimized for large sequences)
    renderSequenceViz(sequence);
    
    // Show results
    document.getElementById('results').classList.remove('hidden');
    
    // Scroll to results
    document.getElementById('results').scrollIntoView({ behavior: 'smooth', block: 'start' });
}

// Load sample sequences
function loadSample(sampleName) {
    const sample = samples[sampleName];
    if (sample) {
        document.getElementById('sequenceInput').value = sample.sequence;
        // Auto-analyze after loading
        setTimeout(() => analyzeSequence(), 300);
    }
}

// Clear all inputs and results
function clearAll() {
    document.getElementById('sequenceInput').value = '';
    document.getElementById('results').classList.add('hidden');
    
    // Remove twin comparison if it exists
    const comparisonDiv = document.getElementById('twinComparison');
    if (comparisonDiv) {
        comparisonDiv.remove();
    }
}

// Parse FASTA format - now handles multiple sequences
function parseFasta(input) {
    const lines = input.trim().split('\n');
    const sequences = [];
    let currentHeader = '';
    let currentSequence = '';
    
    for (let line of lines) {
        line = line.trim();
        if (line.startsWith('>')) {
            // Save previous sequence if it exists
            if (currentHeader && currentSequence) {
                sequences.push({
                    header: currentHeader,
                    sequence: currentSequence
                });
            }
            // Start new sequence
            currentHeader = line.substring(1);
            currentSequence = '';
        } else {
            // Remove any non-DNA characters and convert to uppercase
            currentSequence += line.replace(/[^ATGCatgc]/g, '').toUpperCase();
        }
    }
    
    // Don't forget the last sequence
    if (currentHeader && currentSequence) {
        sequences.push({
            header: currentHeader,
            sequence: currentSequence
        });
    }
    
    // Return single sequence format for backward compatibility
    if (sequences.length === 1) {
        return sequences[0];
    }
    
    // Return multiple sequences
    return sequences.length > 0 ? sequences : { header: '', sequence: '' };
}

// Count nucleotides
function countNucleotides(sequence) {
    const counts = { A: 0, T: 0, G: 0, C: 0 };
    for (let nucleotide of sequence) {
        if (counts.hasOwnProperty(nucleotide)) {
            counts[nucleotide]++;
        }
    }
    return counts;
}

// Calculate GC content
function calculateGC(sequence) {
    const counts = countNucleotides(sequence);
    const total = counts.A + counts.T + counts.G + counts.C;
    const gc = counts.G + counts.C;
    return total > 0 ? (gc / total) * 100 : 0;
}

// Get reverse complement
function reverseComplement(sequence) {
    const complement = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
    };
    
    return sequence
        .split('')
        .reverse()
        .map(base => complement[base] || base)
        .join('');
}

// Find Open Reading Frames (ORFs)
function findORFs(sequence) {
    const orfs = [];
    const startCodon = 'ATG';
    const stopCodons = ['TAA', 'TAG', 'TGA'];
    
    // Search in all three reading frames
    for (let frame = 0; frame < 3; frame++) {
        let i = frame;
        while (i < sequence.length - 2) {
            // Look for start codon
            if (sequence.substr(i, 3) === startCodon) {
                let start = i;
                i += 3;
                
                // Look for stop codon
                while (i < sequence.length - 2) {
                    const codon = sequence.substr(i, 3);
                    if (stopCodons.includes(codon)) {
                        const orfSequence = sequence.substring(start, i + 3);
                        const length = orfSequence.length;
                        
                        // Only include ORFs longer than 30 bp (10 codons)
                        if (length >= 30) {
                            orfs.push({
                                start: start + 1, // 1-indexed
                                end: i + 3,
                                length: length,
                                frame: frame + 1,
                                sequence: orfSequence,
                                protein: translateDNA(orfSequence)
                            });
                        }
                        break;
                    }
                    i += 3;
                }
            } else {
                i += 3;
            }
        }
    }
    
    return orfs;
}

// Translate DNA to protein
function translateDNA(sequence) {
    let protein = '';
    for (let i = 0; i < sequence.length - 2; i += 3) {
        const codon = sequence.substr(i, 3);
        const aa = geneticCode[codon] || 'X';
        if (aa === 'STOP') {
            protein += '*';
            break;
        }
        protein += aa;
    }
    return protein;
}

// Count codons
function countCodons(sequence) {
    const codonCounts = {};
    for (let i = 0; i < sequence.length - 2; i += 3) {
        const codon = sequence.substr(i, 3);
        if (codon.length === 3 && /^[ATGC]{3}$/.test(codon)) {
            codonCounts[codon] = (codonCounts[codon] || 0) + 1;
        }
    }
    return codonCounts;
}

// Render bar chart
function renderBarChart(data, elementId, className = '') {
    const container = document.getElementById(elementId);
    const maxValue = Math.max(...Object.values(data));
    
    let html = '<div class="bar-chart">';
    for (let [label, value] of Object.entries(data)) {
        const height = (value / maxValue) * 100;
        html += `
            <div class="bar-wrapper">
                <div class="bar-value">${value}</div>
                <div class="bar ${className}-${label}" style="height: ${height}%"></div>
                <div class="bar-label">${label}</div>
            </div>
        `;
    }
    html += '</div>';
    
    container.innerHTML = html;
}

// Render GC content chart
function renderGCChart(gcPercent) {
    const container = document.getElementById('gcChart');
    const atPercent = 100 - gcPercent;
    
    // Create a simple conic gradient pie chart
    const pieStyle = `background: conic-gradient(
        var(--nucleotide-g) 0% ${gcPercent/2}%,
        var(--nucleotide-c) ${gcPercent/2}% ${gcPercent}%,
        var(--nucleotide-a) ${gcPercent}% ${gcPercent + atPercent/2}%,
        var(--nucleotide-t) ${gcPercent + atPercent/2}% 100%
    )`;
    
    container.innerHTML = `
        <h3>GC Content Analysis</h3>
        <div class="gc-content-display">
            <div class="pie-chart" style="${pieStyle}"></div>
            <div class="gc-stats">
                <div class="gc-stat-item">
                    <div class="color-box" style="background: linear-gradient(135deg, var(--nucleotide-g), var(--nucleotide-c))"></div>
                    <div>
                        <div class="gc-stat-value">${gcPercent.toFixed(2)}%</div>
                        <div>GC Content</div>
                    </div>
                </div>
                <div class="gc-stat-item">
                    <div class="color-box" style="background: linear-gradient(135deg, var(--nucleotide-a), var(--nucleotide-t))"></div>
                    <div>
                        <div class="gc-stat-value">${atPercent.toFixed(2)}%</div>
                        <div>AT Content</div>
                    </div>
                </div>
            </div>
        </div>
        <div class="explainer mt-20">
            <p><strong>Why is GC content important?</strong> 
            GC-rich regions have stronger hydrogen bonding (3 bonds vs 2 for AT), 
            affecting DNA melting temperature, gene expression, and chromosome structure. 
            Many bacteria have characteristic GC contents (e.g., E. coli ~51%, humans ~41%).</p>
        </div>
    `;
}

// Render sequence visualization
function renderSequenceViz(sequence) {
    const container = document.getElementById('sequenceViz');
    
    // For very large sequences, show summary instead
    if (sequence.length > 100000) {
        container.innerHTML = `
            <h3>Sequence Statistics</h3>
            <div class="stats-grid">
                <div class="stat-item">
                    <span class="stat-value">${(sequence.length/1000000).toFixed(2)}M</span>
                    <span class="stat-label">Total Bases</span>
                </div>
            </div>
            <p style="text-align: center; color: var(--text-secondary); margin-top: 20px;">
                ⚡ Sequence too large for nucleotide-by-nucleotide display<br>
                Use alignment view or download for detailed analysis
            </p>
        `;
        return;
    }
    
    const maxDisplay = 500; // Limit display for performance
    const displaySeq = sequence.substring(0, maxDisplay);
    
    let html = '<h3>Nucleotide-by-Nucleotide View</h3>';
    html += '<div style="line-height: 2;">';
    
    for (let i = 0; i < displaySeq.length; i++) {
        const nucleotide = displaySeq[i];
        html += `<span class="nucleotide ${nucleotide}">${nucleotide}</span>`;
        
        // Add line breaks every 60 characters
        if ((i + 1) % 60 === 0) {
            html += '<br>';
        }
    }
    
    html += '</div>';
    
    if (sequence.length > maxDisplay) {
        html += `<p class="text-center mt-20" style="color: var(--text-secondary);">
            Showing first ${maxDisplay} of ${sequence.length.toLocaleString()} nucleotides
        </p>`;
    }
    
    container.innerHTML = html;
}

// Main analysis function
function analyzeSequence() {
    const input = document.getElementById('sequenceInput').value.trim();
    
    if (!input) {
        alert('Please enter a DNA sequence or load a sample.');
        return;
    }
    
    // Parse FASTA - can return single sequence or array
    const parsed = parseFasta(input);
    
    // Check if we have multiple sequences (twin comparison mode)
    if (Array.isArray(parsed) && parsed.length >= 2) {
        analyzeTwinComparison(parsed);
        return;
    }
    
    // Single sequence mode
    const { header, sequence } = Array.isArray(parsed) ? parsed[0] : parsed;
    
    if (!sequence || sequence.length === 0) {
        alert('No valid DNA sequence found. Please enter A, T, G, C nucleotides.');
        return;
    }
    
    // Perform analyses
    const counts = countNucleotides(sequence);
    const gcContent = calculateGC(sequence);
    const revComp = reverseComplement(sequence);
    const orfs = findORFs(sequence);
    const codonCounts = countCodons(sequence);
    
    // Display basic statistics
    const statsHtml = `
        <div class="stat-item">
            <span class="stat-value">${sequence.length}</span>
            <span class="stat-label">Base Pairs</span>
        </div>
        <div class="stat-item">
            <span class="stat-value" style="color: var(--nucleotide-a)">${counts.A}</span>
            <span class="stat-label">Adenine (A)</span>
        </div>
        <div class="stat-item">
            <span class="stat-value" style="color: var(--nucleotide-t)">${counts.T}</span>
            <span class="stat-label">Thymine (T)</span>
        </div>
        <div class="stat-item">
            <span class="stat-value" style="color: var(--nucleotide-g)">${counts.G}</span>
            <span class="stat-label">Guanine (G)</span>
        </div>
        <div class="stat-item">
            <span class="stat-value" style="color: var(--nucleotide-c)">${counts.C}</span>
            <span class="stat-label">Cytosine (C)</span>
        </div>
        <div class="stat-item">
            <span class="stat-value">${gcContent.toFixed(2)}%</span>
            <span class="stat-label">GC Content</span>
        </div>
    `;
    document.getElementById('basicStats').innerHTML = statsHtml;
    
    // Render charts
    renderBarChart(counts, 'nucleotideChart', 'nucleotide');
    renderGCChart(gcContent);
    
    // Display reverse complement
    document.getElementById('reverseComplement').innerHTML = `
        <h3>Original Sequence (5' → 3'):</h3>
        <div style="margin: 10px 0; word-break: break-all;">${sequence.substring(0, 200)}${sequence.length > 200 ? '...' : ''}</div>
        <h3>Reverse Complement (3' → 5'):</h3>
        <div style="margin: 10px 0; word-break: break-all;">${revComp.substring(0, 200)}${revComp.length > 200 ? '...' : ''}</div>
    `;
    
    // Display ORFs
    if (orfs.length > 0) {
        let orfHtml = '<div class="orf-list">';
        orfs.slice(0, 10).forEach((orf, index) => {
            orfHtml += `
                <div class="orf-item">
                    <div class="orf-header">
                        <div class="orf-title">ORF ${index + 1}</div>
                    </div>
                    <div class="orf-details">
                        <div class="orf-detail"><strong>Position:</strong> ${orf.start} - ${orf.end}</div>
                        <div class="orf-detail"><strong>Length:</strong> ${orf.length} bp (${Math.floor(orf.length/3)} codons)</div>
                        <div class="orf-detail"><strong>Reading Frame:</strong> +${orf.frame}</div>
                    </div>
                    <div class="orf-sequence">
                        <strong>DNA:</strong> ${orf.sequence.substring(0, 90)}${orf.sequence.length > 90 ? '...' : ''}<br>
                        <strong>Protein:</strong> ${orf.protein}
                    </div>
                </div>
            `;
        });
        orfHtml += '</div>';
        
        if (orfs.length > 10) {
            orfHtml += `<p class="text-center mt-20" style="color: var(--text-secondary);">
                Showing 10 of ${orfs.length} ORFs found
            </p>`;
        }
        
        document.getElementById('orfResults').innerHTML = orfHtml;
    } else {
        document.getElementById('orfResults').innerHTML = `
            <p style="color: var(--text-secondary); text-align: center; padding: 20px;">
                No significant ORFs found (minimum length: 30 bp)
            </p>
        `;
    }
    
    // Display codon analysis
    const sortedCodons = Object.entries(codonCounts)
        .sort((a, b) => b[1] - a[1])
        .slice(0, 20);
    
    if (sortedCodons.length > 0) {
        let codonHtml = '<div class="codon-grid">';
        sortedCodons.forEach(([codon, count]) => {
            const aa = geneticCode[codon] || 'Unknown';
            codonHtml += `
                <div class="codon-item">
                    <div class="codon-name">${codon}</div>
                    <div class="codon-count">${aa}</div>
                    <div class="codon-count">Count: ${count}</div>
                </div>
            `;
        });
        codonHtml += '</div>';
        
        if (Object.keys(codonCounts).length > 20) {
            codonHtml += `<p class="text-center mt-20" style="color: var(--text-secondary);">
                Showing top 20 of ${Object.keys(codonCounts).length} codons found
            </p>`;
        }
        
        document.getElementById('codonAnalysis').innerHTML = codonHtml;
    } else {
        document.getElementById('codonAnalysis').innerHTML = `
            <p style="color: var(--text-secondary); text-align: center; padding: 20px;">
                Sequence too short for codon analysis
            </p>
        `;
    }
    
    // Render sequence visualization
    renderSequenceViz(sequence);
    
    // Show results
    document.getElementById('results').classList.remove('hidden');
    
    // Scroll to results
    document.getElementById('results').scrollIntoView({ behavior: 'smooth', block: 'start' });
}

// Twin Comparison Analysis
function analyzeTwinComparison(sequences) {
    // Get first two sequences for comparison
    const twin1 = sequences[0];
    const twin2 = sequences[1];
    
    // Perform analyses on both
    const stats1 = {
        counts: countNucleotides(twin1.sequence),
        gcContent: calculateGC(twin1.sequence),
        length: twin1.sequence.length,
        orfs: findORFs(twin1.sequence),
        codonCounts: countCodons(twin1.sequence)
    };
    
    const stats2 = {
        counts: countNucleotides(twin2.sequence),
        gcContent: calculateGC(twin2.sequence),
        length: twin2.sequence.length,
        orfs: findORFs(twin2.sequence),
        codonCounts: countCodons(twin2.sequence)
    };
    
    // Calculate similarity
    const similarity = calculateSequenceSimilarity(twin1.sequence, twin2.sequence);
    const differences = findSequenceDifferences(twin1.sequence, twin2.sequence);
    
    // Display twin comparison results
    displayTwinComparison(twin1, twin2, stats1, stats2, similarity, differences);
}

// Calculate sequence similarity percentage
function calculateSequenceSimilarity(seq1, seq2) {
    const minLength = Math.min(seq1.length, seq2.length);
    const maxLength = Math.max(seq1.length, seq2.length);
    let matches = 0;
    
    for (let i = 0; i < minLength; i++) {
        if (seq1[i] === seq2[i]) {
            matches++;
        }
    }
    
    return {
        percentage: (matches / maxLength) * 100,
        matches: matches,
        mismatches: minLength - matches,
        lengthDiff: Math.abs(seq1.length - seq2.length)
    };
}

// Find specific differences between sequences and group into variants
function findSequenceDifferences(seq1, seq2) {
    const rawDifferences = [];
    const minLength = Math.min(seq1.length, seq2.length);
    
    // Find all individual differences
    for (let i = 0; i < minLength; i++) {
        if (seq1[i] !== seq2[i]) {
            rawDifferences.push(i);
        }
    }
    
    // Group adjacent differences into variants
    const variants = [];
    let i = 0;
    
    while (i < rawDifferences.length) {
        const startPos = rawDifferences[i];
        let endPos = startPos;
        
        // Check how many consecutive positions differ
        while (i + 1 < rawDifferences.length && rawDifferences[i + 1] === rawDifferences[i] + 1) {
            i++;
            endPos = rawDifferences[i];
        }
        
        const length = endPos - startPos + 1;
        const seq1Segment = seq1.substring(startPos, endPos + 1);
        const seq2Segment = seq2.substring(startPos, endPos + 1);
        
        // Classify variant type
        let variantType;
        let description;
        
        if (length === 1) {
            variantType = 'SNP';
            description = 'Single Nucleotide Polymorphism';
        } else if (length === 2) {
            variantType = 'DNP';
            description = 'Di-Nucleotide Polymorphism';
        } else if (length === 3) {
            variantType = 'TNP';
            description = 'Tri-Nucleotide Polymorphism';
        } else {
            variantType = 'MNP';
            description = `Multi-Nucleotide Polymorphism (${length} bases)`;
        }
        
        // Get larger context
        const contextStart = Math.max(0, startPos - 10);
        const contextEnd = Math.min(seq1.length, endPos + 11);
        const context = seq1.substring(contextStart, contextEnd);
        
        // Highlight the variant in context
        const relativeStart = startPos - contextStart;
        const relativeEnd = endPos - contextStart + 1;
        const contextBefore = context.substring(0, relativeStart);
        const contextVariant = context.substring(relativeStart, relativeEnd);
        const contextAfter = context.substring(relativeEnd);
        
        // Check if this affects a codon boundary (every 3 bases)
        const codonFrame = (startPos % 3) + 1;
        const affectsCodon = Math.floor(startPos / 3) + 1;
        
        // Determine potential impact
        let impact = 'Unknown';
        if (length % 3 === 0) {
            impact = 'In-frame (maintains reading frame)';
        } else {
            impact = 'Frameshift possible';
        }
        
        variants.push({
            type: variantType,
            description: description,
            position: startPos + 1, // 1-indexed
            endPosition: endPos + 1,
            length: length,
            twin1: seq1Segment,
            twin2: seq2Segment,
            context: context,
            contextBefore: contextBefore,
            contextVariant: contextVariant,
            contextAfter: contextAfter,
            codonFrame: codonFrame,
            affectsCodon: affectsCodon,
            impact: impact
        });
        
        i++;
    }
    
    return variants;
}

// Display twin comparison results
function displayTwinComparison(twin1, twin2, stats1, stats2, similarity, differences) {
    // Hide single sequence results, show comparison view
    document.getElementById('results').classList.add('hidden');
    
    // Create comparison view
    let comparisonHTML = `
        <div class="card">
            <h2>👯 Twin Sequence Comparison</h2>
            <div class="explainer">
                <p><strong>Comparing:</strong> ${twin1.header} vs ${twin2.header}</p>
            </div>
            
            <div class="stats-grid">
                <div class="stat-item">
                    <span class="stat-value" style="color: var(--success-color)">${similarity.percentage.toFixed(4)}%</span>
                    <span class="stat-label">Sequence Identity</span>
                </div>
                <div class="stat-item">
                    <span class="stat-value" style="color: var(--primary-color)">${similarity.matches}</span>
                    <span class="stat-label">Matching Bases</span>
                </div>
                <div class="stat-item">
                    <span class="stat-value" style="color: ${differences.length > 0 ? 'var(--warning-color)' : 'var(--success-color)'}">${differences.length}</span>
                    <span class="stat-label">Differences Found</span>
                </div>
                <div class="stat-item">
                    <span class="stat-value" style="color: var(--text-secondary)">${similarity.lengthDiff}</span>
                    <span class="stat-label">Length Difference</span>
                </div>
            </div>
        </div>
        
        <div class="card">
            <h2>📊 Side-by-Side Statistics</h2>
            <div class="twin-comparison-grid">
                <div class="twin-column">
                    <h3 style="color: var(--nucleotide-a)">Twin 1</h3>
                    <p class="twin-header">${twin1.header}</p>
                    <div class="stats-grid">
                        <div class="stat-item">
                            <span class="stat-value">${stats1.length}</span>
                            <span class="stat-label">Base Pairs</span>
                        </div>
                        <div class="stat-item">
                            <span class="stat-value">${stats1.gcContent.toFixed(2)}%</span>
                            <span class="stat-label">GC Content</span>
                        </div>
                        <div class="stat-item">
                            <span class="stat-value">${stats1.orfs.length}</span>
                            <span class="stat-label">ORFs Found</span>
                        </div>
                    </div>
                    <div class="nucleotide-breakdown">
                        <h4>Nucleotide Composition</h4>
                        ${renderNucleotideBreakdown(stats1.counts)}
                    </div>
                </div>
                
                <div class="twin-column">
                    <h3 style="color: var(--nucleotide-t)">Twin 2</h3>
                    <p class="twin-header">${twin2.header}</p>
                    <div class="stats-grid">
                        <div class="stat-item">
                            <span class="stat-value">${stats2.length}</span>
                            <span class="stat-label">Base Pairs</span>
                        </div>
                        <div class="stat-item">
                            <span class="stat-value">${stats2.gcContent.toFixed(2)}%</span>
                            <span class="stat-label">GC Content</span>
                        </div>
                        <div class="stat-item">
                            <span class="stat-value">${stats2.orfs.length}</span>
                            <span class="stat-label">ORFs Found</span>
                        </div>
                    </div>
                    <div class="nucleotide-breakdown">
                        <h4>Nucleotide Composition</h4>
                        ${renderNucleotideBreakdown(stats2.counts)}
                    </div>
                </div>
            </div>
        </div>
    `;
    
    // Add variants/differences section
    if (differences.length > 0) {
        // Count variant types
        const variantCounts = differences.reduce((acc, v) => {
            acc[v.type] = (acc[v.type] || 0) + 1;
            return acc;
        }, {});
        
        const variantSummary = Object.entries(variantCounts)
            .map(([type, count]) => `${count} ${type}${count > 1 ? 's' : ''}`)
            .join(', ');
        
        comparisonHTML += `
            <div class="card">
                <h2>🔍 Genetic Variants Detected</h2>
                <div class="explainer">
                    <p><strong>Found ${differences.length} genetic variant(s):</strong> ${variantSummary}</p>
                    <p>Variants are grouped intelligently - adjacent differences are shown as single events (SNPs, MNPs, etc.)</p>
                </div>
                
                <div class="variant-stats">
                    ${Object.entries(variantCounts).map(([type, count]) => `
                        <div class="variant-stat-item">
                            <div class="variant-count">${count}</div>
                            <div class="variant-label">${type}</div>
                        </div>
                    `).join('')}
                </div>
                
                ${(() => {
                    // Group variants by gene
                    const variantsInGenes = differences.filter(v => findGeneAtPosition(v.position));
                    const variantsByGene = {};
                    
                    variantsInGenes.forEach(variant => {
                        const gene = findGeneAtPosition(variant.position);
                        if (gene) {
                            if (!variantsByGene[gene.name]) {
                                variantsByGene[gene.name] = { gene, variants: [] };
                            }
                            variantsByGene[gene.name].variants.push(variant);
                        }
                    });
                    
                    const geneEntries = Object.entries(variantsByGene).sort((a, b) => b[1].variants.length - a[1].variants.length);
                    
                    if (geneEntries.length === 0) {
                        return '<p class="text-center" style="color: var(--warning-color);">No variants found in known genes</p>';
                    }
                    
                    // Limit display to prevent browser crash
                    const MAX_VARIANTS_PER_GENE = 20;
                    const MAX_GENES_TO_SHOW = 20;
                    
                    let html = `<p class="text-center" style="color: var(--text-secondary); margin-bottom: 20px;">
                        ${variantsInGenes.length} variant${variantsInGenes.length !== 1 ? 's' : ''} across ${geneEntries.length} gene${geneEntries.length !== 1 ? 's' : ''}
                        ${geneEntries.length > MAX_GENES_TO_SHOW ? `<br><span style="color: var(--warning-color);">(Showing top ${MAX_GENES_TO_SHOW} genes)</span>` : ''}
                    </p>`;
                    
                    geneEntries.slice(0, MAX_GENES_TO_SHOW).forEach(([geneName, data]) => {
                        const gene = data.gene;
                        const variants = data.variants.sort((a, b) => b.length - a.length);
                        const variantsToShow = Math.min(MAX_VARIANTS_PER_GENE, variants.length);
                        
                        html += `
                            <div style="margin-bottom: 30px; background: rgba(16, 185, 129, 0.05); padding: 20px; border-radius: 12px; border: 2px solid var(--success-color);">
                                <h4 style="color: var(--success-color); margin-bottom: 10px;">
                                    🧬 ${gene.name} - ${data.variants.length} variant${data.variants.length !== 1 ? 's' : ''}
                                    ${variants.length > MAX_VARIANTS_PER_GENE ? ` <span style="color: var(--warning-color); font-size: 0.9rem;">(showing largest ${MAX_VARIANTS_PER_GENE})</span>` : ''}
                                </h4>
                                <p style="color: var(--text-secondary); font-size: 0.9rem; margin-bottom: 15px;">${gene.description}</p>
                                <div class="differences-list">
                        `;
                        
                        variants.slice(0, MAX_VARIANTS_PER_GENE).forEach((variant, idx) => {
                            const colorClass = variant.type === 'SNP' ? 'snp' : variant.type === 'MNP' ? 'mnp' : 'dnp';
                            html += `
                                <div class="difference-item ${colorClass}">
                                    <div class="diff-header">
                                        <div>
                                            <span class="variant-badge ${variant.type.toLowerCase()}">${variant.type}</span>
                                            <span class="diff-number">Variant ${idx + 1}</span>
                                        </div>
                                        <span class="diff-position">Position: ${variant.position.toLocaleString()}${variant.length > 1 ? `-${variant.endPosition.toLocaleString()}` : ''}</span>
                                    </div>
                                    <div class="diff-content">
                                        <div class="variant-description">
                                            <strong>${variant.description}</strong>
                                        </div>
                                        <div class="diff-bases">
                                            ${renderVariantBases(variant.twin1)}
                                            <span style="margin: 0 15px; font-size: 1.5rem;">→</span>
                                            ${renderVariantBases(variant.twin2)}
                                        </div>
                                        <div class="diff-context">
                                            <strong>Context:</strong> 
                                            <code class="context-display">
                                                ${variant.contextBefore}<span class="variant-highlight">${variant.contextVariant}</span>${variant.contextAfter}
                                            </code>
                                        </div>
                                        <div class="variant-details">
                                            <span>📍 Affects codon #${variant.affectsCodon} (frame ${variant.codonFrame})</span>
                                            <span>⚠️ ${variant.impact}</span>
                                        </div>
                                    </div>
                                </div>
                            `;
                        });
                        
                        html += `
                                </div>
                                ${variants.length > MAX_VARIANTS_PER_GENE ? 
                                    `<div style="text-align: center; padding: 10px; color: var(--text-secondary); font-style: italic;">
                                        ... and ${variants.length - MAX_VARIANTS_PER_GENE} more variant${variants.length - MAX_VARIANTS_PER_GENE !== 1 ? 's' : ''} in this gene
                                    </div>` : ''}
                            </div>
                        `;
                    });
                    
                    // Show hidden genes summary
                    if (geneEntries.length > MAX_GENES_TO_SHOW) {
                        const hiddenGenes = geneEntries.slice(MAX_GENES_TO_SHOW);
                        const hiddenVariantCount = hiddenGenes.reduce((sum, [, data]) => sum + data.variants.length, 0);
                        html += `
                            <div style="text-align: center; padding: 15px; background: rgba(251, 191, 36, 0.1); border-radius: 8px; margin: 20px 0;">
                                <strong>⚠️ ${geneEntries.length - MAX_GENES_TO_SHOW} more gene${geneEntries.length - MAX_GENES_TO_SHOW !== 1 ? 's' : ''} with ${hiddenVariantCount.toLocaleString()} variant${hiddenVariantCount !== 1 ? 's' : ''} not shown</strong>
                            </div>
                        `;
                    }
                    
                    const intergenicCount = differences.length - variantsInGenes.length;
                    if (intergenicCount > 0) {
                        html += `<p class="text-center" style="color: var(--text-secondary); font-style: italic;">
                            (${intergenicCount.toLocaleString()} intergenic variant${intergenicCount !== 1 ? 's' : ''} not shown)
                        </p>`;
                    }
                    
                    return html;
                })()}
            </div>
        `;
    } else {
        comparisonHTML += `
            <div class="card">
                <h2>✅ Perfect Match!</h2>
                <div class="explainer" style="border-color: var(--success-color);">
                    <p>The sequences are 100% identical! This is what we'd expect from identical twins 
                    when looking at the same genomic region.</p>
                </div>
            </div>
        `;
    }
    
    // Add sequence alignment view
    comparisonHTML += `
        <div class="card">
            <h2>🧬 Full Sequence Alignment</h2>
            <div class="explainer">
                <p><strong>Side-by-side view of complete sequences</strong> - Variants are highlighted in color. 
                Scroll horizontally to explore the full length.</p>
            </div>
            <div id="sequenceAlignment"></div>
        </div>
        
        <div class="card">
            <h2>📈 Visual Comparison</h2>
            <div id="comparisonCharts"></div>
        </div>
        
        <div class="card">
            <h2>🧬 Detailed Analysis</h2>
            <div class="button-group">
                <button class="btn-primary" onclick="showDetailedAnalysis(0)">Analyze Twin 1 in Detail</button>
                <button class="btn-primary" onclick="showDetailedAnalysis(1)">Analyze Twin 2 in Detail</button>
            </div>
        </div>
    `;
    
    // Insert comparison view
    const resultsDiv = document.getElementById('results');
    const comparisonDiv = document.createElement('div');
    comparisonDiv.id = 'twinComparison';
    comparisonDiv.innerHTML = comparisonHTML;
    resultsDiv.parentNode.insertBefore(comparisonDiv, resultsDiv);
    
    // Render comparison charts
    renderComparisonCharts(stats1, stats2);
    
    // Render sequence alignment
    renderSequenceAlignment(twin1.sequence, twin2.sequence, differences);
    
    // Store sequences for detailed analysis
    window.currentTwinSequences = [twin1, twin2];
    
    // Scroll to results
    comparisonDiv.scrollIntoView({ behavior: 'smooth', block: 'start' });
}

// Render full sequence alignment view
function renderSequenceAlignment(seq1, seq2, variants) {
    const container = document.getElementById('sequenceAlignment');
    
    // Check if sequence is too large for full rendering
    const isLargeSequence = seq1.length > 100000; // 100kb threshold
    
    if (isLargeSequence) {
        renderLargeSequenceSummary(seq1, seq2, variants, container);
        return;
    }
    
    // Create a map of variant positions for quick lookup
    const variantMap = new Map();
    variants.forEach(variant => {
        for (let i = variant.position - 1; i < variant.endPosition; i++) {
            variantMap.set(i, variant.type);
        }
    });
    
    // Configuration
    const basesPerLine = 60;
    const linesPerBlock = 10;
    const totalLines = Math.ceil(seq1.length / basesPerLine);
    
    // Build alignment HTML
    let html = `
        <div class="alignment-controls">
            <div class="alignment-legend">
                <span class="legend-item"><span class="legend-box match"></span> Match</span>
                <span class="legend-item"><span class="legend-box snp-highlight"></span> SNP</span>
                <span class="legend-item"><span class="legend-box dnp-highlight"></span> DNP</span>
                <span class="legend-item"><span class="legend-box mnp-highlight"></span> MNP/TNP</span>
            </div>
            <div class="alignment-stats">
                <strong>${variants.length}</strong> variant(s) across <strong>${seq1.length}</strong> bases
            </div>
        </div>
        <div class="alignment-container">
    `;
    
    // Generate alignment blocks
    for (let lineStart = 0; lineStart < totalLines; lineStart += linesPerBlock) {
        const lineEnd = Math.min(lineStart + linesPerBlock, totalLines);
        
        html += '<div class="alignment-block">';
        
        for (let line = lineStart; line < lineEnd; line++) {
            const start = line * basesPerLine;
            const end = Math.min(start + basesPerLine, seq1.length);
            const position = start + 1; // 1-indexed
            
            // Position ruler
            html += `<div class="alignment-line">`;
            html += `<div class="position-label">${position}</div>`;
            html += `<div class="sequence-row">`;
            
            // Twin 1 sequence
            for (let i = start; i < end; i++) {
                const base = seq1[i];
                const variantType = variantMap.get(i);
                const isMatch = seq1[i] === seq2[i];
                
                let className = 'align-base';
                if (variantType) {
                    className += ` variant-${variantType.toLowerCase()}`;
                } else if (isMatch) {
                    className += ' match';
                }
                
                html += `<span class="${className}" data-pos="${i + 1}">${base}</span>`;
            }
            
            html += `</div></div>`;
            
            // Match indicator line
            html += `<div class="alignment-line">`;
            html += `<div class="position-label"></div>`;
            html += `<div class="sequence-row match-indicator">`;
            
            for (let i = start; i < end; i++) {
                const isMatch = seq1[i] === seq2[i];
                const symbol = isMatch ? '|' : '×';
                const className = isMatch ? 'match-symbol' : 'diff-symbol';
                html += `<span class="${className}">${symbol}</span>`;
            }
            
            html += `</div></div>`;
            
            // Twin 2 sequence
            html += `<div class="alignment-line">`;
            html += `<div class="position-label">${position}</div>`;
            html += `<div class="sequence-row">`;
            
            for (let i = start; i < end; i++) {
                const base = seq2[i];
                const variantType = variantMap.get(i);
                const isMatch = seq1[i] === seq2[i];
                
                let className = 'align-base';
                if (variantType) {
                    className += ` variant-${variantType.toLowerCase()}`;
                } else if (isMatch) {
                    className += ' match';
                }
                
                html += `<span class="${className}" data-pos="${i + 1}">${base}</span>`;
            }
            
            html += `</div></div>`;
            
            // Spacer between groups
            if (line < lineEnd - 1) {
                html += '<div class="alignment-spacer"></div>';
            }
        }
        
        html += '</div>'; // alignment-block
    }
    
    html += '</div>'; // alignment-container
    
    // Add variant density map
    html += renderVariantDensityMap(seq1.length, variants);
    
    container.innerHTML = html;
    
    // Add click handlers for highlighting
    addAlignmentInteractivity();
}

// Render summary view for large sequences
function renderLargeSequenceSummary(seq1, seq2, variants, container) {
    const seqLength = seq1.length;
    const diffs = variants.reduce((sum, v) => sum + v.length, 0);
    const similarity = ((seqLength - diffs) / seqLength * 100).toFixed(4);
    
    let html = `
        <div class="large-sequence-notice">
            <h3>⚡ Large Sequence Mode</h3>
            <p><strong>Sequence too large for full display (${(seqLength/1000000).toFixed(2)}M bases)</strong></p>
            <p>Showing summary statistics and variant density map instead.</p>
        </div>
        
        <div class="stats-grid" style="margin: 20px 0;">
            <div class="stat-item">
                <span class="stat-value">${(seqLength/1000000).toFixed(2)}M</span>
                <span class="stat-label">Total Bases</span>
            </div>
            <div class="stat-item">
                <span class="stat-value">${diffs.toLocaleString()}</span>
                <span class="stat-label">Differences</span>
            </div>
            <div class="stat-item">
                <span class="stat-value">${similarity}%</span>
                <span class="stat-label">Similarity</span>
            </div>
            <div class="stat-item">
                <span class="stat-value">${variants.length}</span>
                <span class="stat-label">Variants</span>
            </div>
        </div>
    `;
    
    // Add variant density map
    html += renderVariantDensityMap(seqLength, variants);
    
    // Add sample regions around variants
    html += `
        <div style="margin-top: 30px;">
            <h4>Largest Variants (Top 100 by Size)</h4>
            <p style="color: var(--text-secondary); margin-bottom: 15px;">
                Showing the most significant genetic differences with gene annotations. 
                <a href="#" onclick="downloadFullAlignment(); return false;">Download full data</a> for complete analysis.
            </p>
        </div>
    `;
    
    // Filter to only show variants in known genes
    const variantsInGenes = variants.filter(v => findGeneAtPosition(v.position));
    
    if (variantsInGenes.length === 0) {
        html += `
            <div style="background: rgba(251, 191, 36, 0.1); padding: 15px; border-radius: 8px; margin-bottom: 20px; border-left: 4px solid var(--warning-color);">
                <strong>⚠️ No variants found in known genes</strong><br>
                <span style="color: var(--text-secondary); font-size: 0.9rem;">
                    This could mean variants are in intergenic regions, or position information needs adjustment.
                    Total variants found: ${variants.length.toLocaleString()}
                </span>
            </div>
        `;
        
        // Show first 10 variants anyway for debugging
        html += `<h4>Showing first 10 variants (for reference):</h4>`;
        sortedVariants.slice(0, 10).forEach((variant, idx) => {
            const pos = variant.position - 1;
            const contextStart = Math.max(0, pos - 30);
            const contextEnd = Math.min(seq1.length, pos + 30);
            
            const context1 = seq1.substring(contextStart, contextEnd);
            const context2 = seq2.substring(contextStart, contextEnd);
            
            html += `
                <div class="sample-region">
                    <div style="font-weight: bold; margin-bottom: 8px;">
                        #${idx + 1} ${variant.type} at position ${variant.position.toLocaleString()} (${variant.length} bases)
                    </div>
                    <div style="font-family: 'Courier New', monospace; font-size: 0.85rem; line-height: 1.6;">
                        <div>Twin 1: ${highlightVariantInContext(context1, pos - contextStart, variant.length)}</div>
                        <div>Twin 2: ${highlightVariantInContext(context2, pos - contextStart, variant.length)}</div>
                    </div>
                </div>
            `;
        });
        
        container.innerHTML = html;
        return;
    }
    
    // Group variants by gene
    const variantsByGene = {};
    variantsInGenes.forEach(variant => {
        const gene = findGeneAtPosition(variant.position);
        if (gene) {
            if (!variantsByGene[gene.name]) {
                variantsByGene[gene.name] = {
                    gene: gene,
                    variants: []
                };
            }
            variantsByGene[gene.name].variants.push(variant);
        }
    });
    
    // Sort genes by number of variants (most variants first)
    const geneEntries = Object.entries(variantsByGene).sort((a, b) => b[1].variants.length - a[1].variants.length);
    
    // Limit display to prevent browser crash with huge datasets
    const MAX_VARIANTS_PER_GENE = 20;
    const MAX_GENES_TO_SHOW = 20;
    
    // Show summary
    const genesToShow = Math.min(MAX_GENES_TO_SHOW, geneEntries.length);
    html += `
        <div style="background: rgba(16, 185, 129, 0.1); padding: 15px; border-radius: 8px; margin-bottom: 20px; border-left: 4px solid var(--success-color);">
            <strong>🧬 Variants Grouped by Gene:</strong><br>
            Found ${variantsInGenes.length.toLocaleString()} variant${variantsInGenes.length !== 1 ? 's' : ''} across ${geneEntries.length} gene${geneEntries.length !== 1 ? 's' : ''}
            ${geneEntries.length > MAX_GENES_TO_SHOW ? `<br><span style="color: var(--text-secondary); font-size: 0.9rem;">(Showing top ${MAX_GENES_TO_SHOW} genes with most variants)</span>` : ''}
        </div>
    `;
    
    // Show gene summary table (limit to prevent crash)
    html += `
        <div style="background: rgba(0, 0, 0, 0.3); padding: 20px; border-radius: 8px; margin-bottom: 25px;">
            <h4 style="margin-bottom: 15px;">Gene Summary</h4>
            <div style="display: grid; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); gap: 15px;">
                ${geneEntries.slice(0, MAX_GENES_TO_SHOW).map(([geneName, data]) => `
                    <div style="background: rgba(16, 185, 129, 0.1); padding: 12px; border-radius: 6px; border-left: 3px solid var(--success-color);">
                        <div style="font-weight: bold; color: var(--success-color);">${geneName}</div>
                        <div style="font-size: 0.85rem; color: var(--text-secondary); margin-top: 4px;">
                            ${data.variants.length} variant${data.variants.length !== 1 ? 's' : ''}
                        </div>
                    </div>
                `).join('')}
            </div>
        </div>
    `;
    
    // Display variants grouped by gene (with limits)
    geneEntries.slice(0, MAX_GENES_TO_SHOW).forEach(([geneName, data]) => {
        const gene = data.gene;
        const geneVariants = data.variants.sort((a, b) => b.length - a.length); // Sort by size within gene
        const variantsToShow = Math.min(MAX_VARIANTS_PER_GENE, geneVariants.length);
        
        html += `
            <div style="background: linear-gradient(135deg, rgba(16, 185, 129, 0.15), rgba(16, 185, 129, 0.05)); padding: 20px; border-radius: 12px; margin-bottom: 25px; border: 2px solid var(--success-color);">
                <div style="margin-bottom: 20px;">
                    <h3 style="color: var(--success-color); margin-bottom: 8px;">
                        🧬 ${gene.name} Gene
                    </h3>
                    <div style="color: var(--text-secondary); margin-bottom: 8px;">
                        ${gene.description}
                    </div>
                    <div style="display: flex; gap: 20px; flex-wrap: wrap; font-size: 0.9rem; color: var(--text-secondary);">
                        <span>📍 Position: ${gene.start.toLocaleString()}-${gene.end.toLocaleString()}</span>
                        <span>📏 Length: ${((gene.end - gene.start) / 1000).toFixed(1)}kb</span>
                        <span>🔬 Variants: ${geneVariants.length}</span>
                        ${geneVariants.length > MAX_VARIANTS_PER_GENE ? `<span style="color: var(--warning-color);">⚠️ Showing largest ${MAX_VARIANTS_PER_GENE}</span>` : ''}
                    </div>
                </div>
                
                <div style="display: flex; flex-direction: column; gap: 15px;">
        `;
        
        geneVariants.slice(0, MAX_VARIANTS_PER_GENE).forEach((variant, idx) => {
            const pos = variant.position - 1;
            const contextStart = Math.max(0, pos - 30);
            const contextEnd = Math.min(seq1.length, pos + 30);
            
            const context1 = seq1.substring(contextStart, contextEnd);
            const context2 = seq2.substring(contextStart, contextEnd);
            
            const sizeLabel = variant.length >= 4 ? '⚠️' : variant.length === 3 ? '🔸' : variant.length === 2 ? '🔹' : '•';
            
            html += `
                <div style="background: rgba(0, 0, 0, 0.3); padding: 15px; border-radius: 8px; border-left: 3px solid ${
                    variant.type === 'MNP' ? 'var(--danger-color)' : 
                    variant.type === 'TNP' ? 'var(--warning-color)' : 
                    variant.type === 'DNP' ? '#f97316' : 'var(--primary-color)'
                };">
                    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px;">
                        <div style="font-weight: bold;">
                            ${sizeLabel} ${variant.type} at position ${variant.position.toLocaleString()}
                        </div>
                        <div style="color: var(--text-secondary); font-size: 0.9rem;">
                            ${variant.length} base${variant.length > 1 ? 's' : ''}
                        </div>
                    </div>
                    <div style="font-family: 'Courier New', monospace; font-size: 0.85rem; line-height: 1.6;">
                        <div>Twin 1: ${highlightVariantInContext(context1, pos - contextStart, variant.length)}</div>
                        <div>Twin 2: ${highlightVariantInContext(context2, pos - contextStart, variant.length)}</div>
                    </div>
                </div>
            `;
        });
        
        html += `
                </div>
                ${geneVariants.length > MAX_VARIANTS_PER_GENE ? 
                    `<div style="text-align: center; padding: 15px; color: var(--text-secondary); font-style: italic;">
                        ... and ${geneVariants.length - MAX_VARIANTS_PER_GENE} more variant${geneVariants.length - MAX_VARIANTS_PER_GENE !== 1 ? 's' : ''} in this gene
                    </div>` : ''}
            </div>
        `;
    });
    
    // Show summary of hidden genes
    if (geneEntries.length > MAX_GENES_TO_SHOW) {
        const hiddenGenes = geneEntries.slice(MAX_GENES_TO_SHOW);
        const hiddenVariantCount = hiddenGenes.reduce((sum, [, data]) => sum + data.variants.length, 0);
        html += `
            <div style="text-align: center; padding: 20px; background: rgba(251, 191, 36, 0.1); border-radius: 8px; border-left: 4px solid var(--warning-color); margin-top: 20px;">
                <strong>⚠️ ${geneEntries.length - MAX_GENES_TO_SHOW} more gene${geneEntries.length - MAX_GENES_TO_SHOW !== 1 ? 's' : ''} with ${hiddenVariantCount.toLocaleString()} variant${hiddenVariantCount !== 1 ? 's' : ''} not shown</strong><br>
                <span style="color: var(--text-secondary); font-size: 0.9rem;">
                    (Limiting display to prevent browser performance issues. Download full data for complete analysis.)
                </span>
            </div>
        `;
    }
    
    // Show total intergenic variants
    const intergenicCount = variants.length - variantsInGenes.length;
    if (intergenicCount > 0) {
        html += `<div style="text-align: center; color: var(--text-secondary); margin-top: 20px; padding: 15px; background: rgba(100, 100, 100, 0.1); border-radius: 8px;">
            <span style="font-style: italic;">
                ${intergenicCount.toLocaleString()} additional variant${intergenicCount !== 1 ? 's' : ''} in intergenic regions (not shown)
            </span>
        </div>`;
    }
    
    container.innerHTML = html;
}

// Helper to highlight variant in context
function highlightVariantInContext(context, variantOffset, variantLength) {
    const before = context.substring(0, variantOffset);
    const variant = context.substring(variantOffset, variantOffset + variantLength);
    const after = context.substring(variantOffset + variantLength);
    
    return `${before}<span style="background: var(--warning-color); color: #000; padding: 2px 4px; border-radius: 3px; font-weight: bold;">${variant}</span>${after}`;
}

// Download full alignment data
function downloadFullAlignment() {
    if (!window.currentTwinSequences || window.currentTwinSequences.length < 2) {
        alert('No twin data available');
        return;
    }
    
    const seq1 = window.currentTwinSequences[0];
    const seq2 = window.currentTwinSequences[1];
    
    const content = `>${seq1.header}\n${seq1.sequence}\n\n>${seq2.header}\n${seq2.sequence}`;
    
    const blob = new Blob([content], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'twin_sequences_full.fasta';
    a.click();
    URL.revokeObjectURL(url);
    
    alert('Full sequences downloaded as FASTA file!');
}

// Render a visual map showing variant density across the sequence
function renderVariantDensityMap(seqLength, variants) {
    const buckets = 100; // Divide sequence into 100 segments
    const bucketSize = Math.ceil(seqLength / buckets);
    const densityMap = new Array(buckets).fill(0);
    
    // Count variants in each bucket
    variants.forEach(variant => {
        const bucket = Math.floor((variant.position - 1) / bucketSize);
        if (bucket < buckets) {
            densityMap[bucket]++;
        }
    });
    
    const maxDensity = Math.max(...densityMap, 1);
    
    let html = `
        <div class="density-map-container">
            <h4>Variant Density Map</h4>
            <div class="density-map">
    `;
    
    for (let i = 0; i < buckets; i++) {
        const density = densityMap[i];
        const height = (density / maxDensity) * 100;
        const position = i * bucketSize + 1;
        const hasVariant = density > 0;
        
        html += `
            <div class="density-bar ${hasVariant ? 'has-variant' : ''}" 
                 style="height: ${height}%"
                 title="Position ${position}-${Math.min(position + bucketSize - 1, seqLength)}: ${density} variant(s)"
                 data-bucket="${i}">
            </div>
        `;
    }
    
    html += `
            </div>
            <div class="density-labels">
                <span>Start</span>
                <span>Position</span>
                <span>End</span>
            </div>
        </div>
    `;
    
    return html;
}

// Add interactivity to alignment view
function addAlignmentInteractivity() {
    // Hover effect to highlight corresponding positions
    const bases = document.querySelectorAll('.align-base');
    bases.forEach(base => {
        base.addEventListener('mouseenter', function() {
            const pos = this.getAttribute('data-pos');
            document.querySelectorAll(`[data-pos="${pos}"]`).forEach(b => {
                b.classList.add('highlight-hover');
            });
        });
        
        base.addEventListener('mouseleave', function() {
            document.querySelectorAll('.highlight-hover').forEach(b => {
                b.classList.remove('highlight-hover');
            });
        });
        
        // Click to show position info
        base.addEventListener('click', function() {
            const pos = this.getAttribute('data-pos');
            showPositionInfo(pos);
        });
    });
    
    // Make density map interactive
    const densityBars = document.querySelectorAll('.density-bar');
    densityBars.forEach(bar => {
        bar.addEventListener('click', function() {
            const bucket = parseInt(this.getAttribute('data-bucket'));
            scrollToAlignmentBucket(bucket);
        });
    });
}

// Show information about a specific position
function showPositionInfo(position) {
    const pos = parseInt(position) - 1; // Convert to 0-indexed
    
    if (!window.currentTwinSequences || window.currentTwinSequences.length < 2) {
        return;
    }
    
    const seq1 = window.currentTwinSequences[0].sequence;
    const seq2 = window.currentTwinSequences[1].sequence;
    
    const base1 = seq1[pos];
    const base2 = seq2[pos];
    const codon = Math.floor(pos / 3) + 1;
    const frame = (pos % 3) + 1;
    const isMatch = base1 === base2;
    
    const info = `
        <strong>Position ${position}</strong><br>
        Twin 1: <span class="nucleotide ${base1}">${base1}</span><br>
        Twin 2: <span class="nucleotide ${base2}">${base2}</span><br>
        ${isMatch ? '✅ Match' : '⚠️ Difference'}<br>
        Codon: #${codon}, Frame: ${frame}
    `;
    
    // Create or update tooltip
    let tooltip = document.getElementById('alignment-tooltip');
    if (!tooltip) {
        tooltip = document.createElement('div');
        tooltip.id = 'alignment-tooltip';
        tooltip.className = 'alignment-tooltip';
        document.body.appendChild(tooltip);
    }
    
    tooltip.innerHTML = info;
    tooltip.style.display = 'block';
    
    // Position near cursor
    const clickedElement = document.querySelector(`[data-pos="${position}"]`);
    if (clickedElement) {
        const rect = clickedElement.getBoundingClientRect();
        tooltip.style.left = (rect.left + window.scrollX) + 'px';
        tooltip.style.top = (rect.bottom + window.scrollY + 5) + 'px';
    }
    
    // Hide after 3 seconds
    setTimeout(() => {
        tooltip.style.display = 'none';
    }, 3000);
}

// Scroll to a specific bucket in the alignment
function scrollToAlignmentBucket(bucket) {
    const blocks = document.querySelectorAll('.alignment-block');
    const targetBlock = Math.floor(bucket / 10);
    
    if (blocks[targetBlock]) {
        blocks[targetBlock].scrollIntoView({ behavior: 'smooth', block: 'center' });
    }
}

// Helper function to render variant bases with proper spacing
function renderVariantBases(bases) {
    if (bases.length === 1) {
        return `<span class="nucleotide ${bases}">${bases}</span>`;
    }
    // Multiple bases - show each one
    return bases.split('').map(base => 
        `<span class="nucleotide ${base}" style="margin: 0 2px;">${base}</span>`
    ).join('');
}

// Helper function to render nucleotide breakdown
function renderNucleotideBreakdown(counts) {
    return `
        <div style="display: flex; justify-content: space-around; margin-top: 10px;">
            <div style="text-align: center;">
                <div class="nucleotide A" style="padding: 10px; font-size: 1.2rem;">${counts.A}</div>
                <div style="font-size: 0.8rem; margin-top: 5px;">A</div>
            </div>
            <div style="text-align: center;">
                <div class="nucleotide T" style="padding: 10px; font-size: 1.2rem;">${counts.T}</div>
                <div style="font-size: 0.8rem; margin-top: 5px;">T</div>
            </div>
            <div style="text-align: center;">
                <div class="nucleotide G" style="padding: 10px; font-size: 1.2rem;">${counts.G}</div>
                <div style="font-size: 0.8rem; margin-top: 5px;">G</div>
            </div>
            <div style="text-align: center;">
                <div class="nucleotide C" style="padding: 10px; font-size: 1.2rem;">${counts.C}</div>
                <div style="font-size: 0.8rem; margin-top: 5px;">C</div>
            </div>
        </div>
    `;
}

// Render comparison charts
function renderComparisonCharts(stats1, stats2) {
    const container = document.getElementById('comparisonCharts');
    
    let html = '<h3>Nucleotide Composition Comparison</h3>';
    html += '<div style="display: grid; grid-template-columns: 1fr 1fr; gap: 20px;">';
    
    // Twin 1 chart
    html += '<div><h4 style="color: var(--nucleotide-a); text-align: center;">Twin 1</h4>';
    html += '<div id="chart1"></div></div>';
    
    // Twin 2 chart
    html += '<div><h4 style="color: var(--nucleotide-t); text-align: center;">Twin 2</h4>';
    html += '<div id="chart2"></div></div>';
    
    html += '</div>';
    
    // GC Content comparison
    html += '<h3 style="margin-top: 30px;">GC Content Comparison</h3>';
    html += '<div style="display: grid; grid-template-columns: 1fr 1fr; gap: 20px;">';
    html += `
        <div style="text-align: center;">
            <div style="font-size: 3rem; font-weight: bold; color: var(--nucleotide-g);">${stats1.gcContent.toFixed(2)}%</div>
            <div>Twin 1 GC Content</div>
        </div>
        <div style="text-align: center;">
            <div style="font-size: 3rem; font-weight: bold; color: var(--nucleotide-c);">${stats2.gcContent.toFixed(2)}%</div>
            <div>Twin 2 GC Content</div>
        </div>
    `;
    html += '</div>';
    
    container.innerHTML = html;
    
    // Render individual bar charts
    renderBarChart(stats1.counts, 'chart1', 'nucleotide');
    renderBarChart(stats2.counts, 'chart2', 'nucleotide');
}

// Show detailed analysis for a specific twin
function showDetailedAnalysis(twinIndex) {
    if (!window.currentTwinSequences || !window.currentTwinSequences[twinIndex]) {
        alert('Twin sequence data not available');
        return;
    }
    
    const twin = window.currentTwinSequences[twinIndex];
    
    // Hide comparison view
    const comparisonDiv = document.getElementById('twinComparison');
    if (comparisonDiv) {
        comparisonDiv.style.display = 'none';
    }
    
    // Create single sequence input format
    const fastaInput = `>${twin.header}\n${twin.sequence}`;
    document.getElementById('sequenceInput').value = fastaInput;
    
    // Run single sequence analysis
    const { header, sequence } = twin;
    
    if (!sequence || sequence.length === 0) {
        alert('No valid DNA sequence found.');
        return;
    }
    
    // Perform analyses
    const counts = countNucleotides(sequence);
    const gcContent = calculateGC(sequence);
    const revComp = reverseComplement(sequence);
    const orfs = findORFs(sequence);
    const codonCounts = countCodons(sequence);
    
    // Display basic statistics
    const statsHtml = `
        <div class="stat-item">
            <span class="stat-value">${sequence.length}</span>
            <span class="stat-label">Base Pairs</span>
        </div>
        <div class="stat-item">
            <span class="stat-value" style="color: var(--nucleotide-a)">${counts.A}</span>
            <span class="stat-label">Adenine (A)</span>
        </div>
        <div class="stat-item">
            <span class="stat-value" style="color: var(--nucleotide-t)">${counts.T}</span>
            <span class="stat-label">Thymine (T)</span>
        </div>
        <div class="stat-item">
            <span class="stat-value" style="color: var(--nucleotide-g)">${counts.G}</span>
            <span class="stat-label">Guanine (G)</span>
        </div>
        <div class="stat-item">
            <span class="stat-value" style="color: var(--nucleotide-c)">${counts.C}</span>
            <span class="stat-label">Cytosine (C)</span>
        </div>
        <div class="stat-item">
            <span class="stat-value">${gcContent.toFixed(2)}%</span>
            <span class="stat-label">GC Content</span>
        </div>
    `;
    document.getElementById('basicStats').innerHTML = statsHtml;
    
    // Render charts
    renderBarChart(counts, 'nucleotideChart', 'nucleotide');
    renderGCChart(gcContent);
    
    // Display reverse complement
    document.getElementById('reverseComplement').innerHTML = `
        <h3>Original Sequence (5' → 3'):</h3>
        <div style="margin: 10px 0; word-break: break-all;">${sequence.substring(0, 200)}${sequence.length > 200 ? '...' : ''}</div>
        <h3>Reverse Complement (3' → 5'):</h3>
        <div style="margin: 10px 0; word-break: break-all;">${revComp.substring(0, 200)}${revComp.length > 200 ? '...' : ''}</div>
    `;
    
    // Display ORFs
    if (orfs.length > 0) {
        let orfHtml = '<div class="orf-list">';
        orfs.slice(0, 10).forEach((orf, index) => {
            orfHtml += `
                <div class="orf-item">
                    <div class="orf-header">
                        <div class="orf-title">ORF ${index + 1}</div>
                    </div>
                    <div class="orf-details">
                        <div class="orf-detail"><strong>Position:</strong> ${orf.start} - ${orf.end}</div>
                        <div class="orf-detail"><strong>Length:</strong> ${orf.length} bp (${Math.floor(orf.length/3)} codons)</div>
                        <div class="orf-detail"><strong>Reading Frame:</strong> +${orf.frame}</div>
                    </div>
                    <div class="orf-sequence">
                        <strong>DNA:</strong> ${orf.sequence.substring(0, 90)}${orf.sequence.length > 90 ? '...' : ''}<br>
                        <strong>Protein:</strong> ${orf.protein}
                    </div>
                </div>
            `;
        });
        orfHtml += '</div>';
        
        if (orfs.length > 10) {
            orfHtml += `<p class="text-center mt-20" style="color: var(--text-secondary);">
                Showing 10 of ${orfs.length} ORFs found
            </p>`;
        }
        
        document.getElementById('orfResults').innerHTML = orfHtml;
    } else {
        document.getElementById('orfResults').innerHTML = `
            <p style="color: var(--text-secondary); text-align: center; padding: 20px;">
                No significant ORFs found (minimum length: 30 bp)
            </p>
        `;
    }
    
    // Display codon analysis
    const sortedCodons = Object.entries(codonCounts)
        .sort((a, b) => b[1] - a[1])
        .slice(0, 20);
    
    if (sortedCodons.length > 0) {
        let codonHtml = '<div class="codon-grid">';
        sortedCodons.forEach(([codon, count]) => {
            const aa = geneticCode[codon] || 'Unknown';
            codonHtml += `
                <div class="codon-item">
                    <div class="codon-name">${codon}</div>
                    <div class="codon-count">${aa}</div>
                    <div class="codon-count">Count: ${count}</div>
                </div>
            `;
        });
        codonHtml += '</div>';
        
        if (Object.keys(codonCounts).length > 20) {
            codonHtml += `<p class="text-center mt-20" style="color: var(--text-secondary);">
                Showing top 20 of ${Object.keys(codonCounts).length} codons found
            </p>`;
        }
        
        document.getElementById('codonAnalysis').innerHTML = codonHtml;
    } else {
        document.getElementById('codonAnalysis').innerHTML = `
            <p style="color: var(--text-secondary); text-align: center; padding: 20px;">
                Sequence too short for codon analysis
            </p>
        `;
    }
    
    // Render sequence visualization
    renderSequenceViz(sequence);
    
    // Show results
    document.getElementById('results').classList.remove('hidden');
    
    // Scroll to results
    document.getElementById('results').scrollIntoView({ behavior: 'smooth', block: 'start' });
}


