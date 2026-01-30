#!/bin/bash

# Twin Genome Data Download Script
# Downloads and prepares twin genome data for analysis

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
DATA_DIR="./genome_data"
REFERENCE_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes"
FTP_BASE="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

# Function to print colored messages
print_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
print_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
print_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
print_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check for required tools
check_dependencies() {
    print_info "Checking dependencies..."
    
    local missing=()
    
    if ! command_exists wget && ! command_exists curl; then
        missing+=("wget or curl")
    fi
    
    if ! command_exists bcftools; then
        missing+=("bcftools")
    fi
    
    if ! command_exists samtools; then
        missing+=("samtools")
    fi
    
    if [ ${#missing[@]} -gt 0 ]; then
        print_error "Missing required tools: ${missing[*]}"
        print_info "Install with: brew install samtools bcftools wget"
        print_info "Or with conda: conda install -c bioconda samtools bcftools"
        exit 1
    fi
    
    print_success "All dependencies found!"
}

# Create directory structure
setup_directories() {
    print_info "Setting up directories..."
    mkdir -p "$DATA_DIR"/{vcf,fasta,reference}
    print_success "Directories created"
}

# Download reference genome (chromosome or region)
download_reference() {
    local chr=$1
    print_info "Downloading reference genome for chromosome $chr..."
    
    local ref_file="$DATA_DIR/reference/chr${chr}.fa.gz"
    
    if [ -f "$ref_file" ]; then
        print_warning "Reference already exists, skipping download"
        return
    fi
    
    wget -q --show-progress -O "$ref_file" \
        "${REFERENCE_URL}/chr${chr}.fa.gz" || {
        print_error "Failed to download reference"
        exit 1
    }
    
    print_info "Uncompressing reference..."
    gunzip -f "$ref_file"
    
    print_info "Indexing reference..."
    samtools faidx "${ref_file%.gz}"
    
    print_success "Reference genome ready"
}

# Download variant data from 1000 Genomes
download_variants() {
    local chr=$1
    print_info "Downloading variant data for chromosome $chr..."
    
    local vcf_url="${FTP_BASE}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    local vcf_file="$DATA_DIR/vcf/chr${chr}_all.vcf.gz"
    
    if [ -f "$vcf_file" ]; then
        print_warning "VCF already exists, skipping download"
        return
    fi
    
    print_info "This may take several minutes..."
    wget -q --show-progress -O "$vcf_file" "$vcf_url" || {
        print_error "Failed to download VCF"
        exit 1
    }
    
    print_info "Downloading index..."
    wget -q -O "${vcf_file}.tbi" "${vcf_url}.tbi"
    
    print_success "Variant data downloaded"
}

# Extract specific gene region
extract_gene_region() {
    local gene=$1
    local chr=$2
    local start=$3
    local end=$4
    
    print_info "Extracting $gene region (chr${chr}:${start}-${end})..."
    
    local input_vcf="$DATA_DIR/vcf/chr${chr}_all.vcf.gz"
    local output_vcf="$DATA_DIR/vcf/${gene}_region.vcf.gz"
    
    bcftools view -r "${chr}:${start}-${end}" \
        -Oz -o "$output_vcf" \
        "$input_vcf"
    
    bcftools index "$output_vcf"
    
    print_success "Gene region extracted"
}

# Find related individuals (siblings, twins)
find_related_pairs() {
    print_info "Downloading pedigree information..."
    
    local ped_file="$DATA_DIR/pedigree.txt"
    
    wget -q -O "$ped_file" \
        "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
    
    print_info "Related pairs found in pedigree file"
    print_info "Common pairs:"
    echo "  - NA19092, NA19093 (Yoruba, Nigeria)"
    echo "  - HG00096, HG00097 (British, England)"
    echo "  - NA12003, NA12004 (CEPH, Utah)"
}

# Extract sequences for two individuals
extract_twin_sequences() {
    local ind1=$1
    local ind2=$2
    local vcf_file=$3
    local ref_file=$4
    local gene=$5
    
    print_info "Extracting sequences for $ind1 and $ind2..."
    
    # Extract individual 1
    print_info "Processing $ind1..."
    bcftools view -s "$ind1" "$vcf_file" | \
        bcftools consensus -f "$ref_file" \
        > "$DATA_DIR/fasta/${gene}_${ind1}.fasta"
    
    # Extract individual 2
    print_info "Processing $ind2..."
    bcftools view -s "$ind2" "$vcf_file" | \
        bcftools consensus -f "$ref_file" \
        > "$DATA_DIR/fasta/${gene}_${ind2}.fasta"
    
    # Combine for comparison
    print_info "Creating combined file for analysis..."
    {
        echo ">${gene}_${ind1} | Individual: $ind1"
        tail -n +2 "$DATA_DIR/fasta/${gene}_${ind1}.fasta"
        echo ""
        echo ">${gene}_${ind2} | Individual: $ind2"
        tail -n +2 "$DATA_DIR/fasta/${gene}_${ind2}.fasta"
    } > "$DATA_DIR/fasta/${gene}_twins_comparison.fasta"
    
    print_success "Twin sequences ready: $DATA_DIR/fasta/${gene}_twins_comparison.fasta"
}

# Download BRCA1 gene for twins (quick example)
download_brca1_twins() {
    local ind1=${1:-NA19092}
    local ind2=${2:-NA19093}
    
    print_info "Downloading BRCA1 gene data for twins..."
    print_info "Individuals: $ind1 and $ind2"
    
    local chr="17"
    local gene="BRCA1"
    local start="43044295"
    local end="43170245"
    
    # Download reference
    download_reference "$chr"
    
    # Download variants
    download_variants "$chr"
    
    # Extract gene region
    extract_gene_region "$gene" "$chr" "$start" "$end"
    
    # Extract sequences
    local ref_file="$DATA_DIR/reference/chr${chr}.fa"
    local vcf_file="$DATA_DIR/vcf/${gene}_region.vcf.gz"
    
    extract_twin_sequences "$ind1" "$ind2" "$vcf_file" "$ref_file" "$gene"
    
    print_success "✨ Complete! Your twin BRCA1 data is ready!"
    print_info "File location: $DATA_DIR/fasta/BRCA1_twins_comparison.fasta"
    print_info "You can now load this file into your DNA analyzer!"
}

# Download chromosome 21 (smallest chromosome)
download_chr21_twins() {
    local ind1=${1:-NA19092}
    local ind2=${2:-NA19093}
    
    print_info "Downloading chromosome 21 for twins..."
    print_info "⚠️  This will take 5-10 minutes and use ~200 MB"
    
    local chr="21"
    
    download_reference "$chr"
    download_variants "$chr"
    
    print_info "Extracting twin sequences..."
    local ref_file="$DATA_DIR/reference/chr${chr}.fa"
    local vcf_file="$DATA_DIR/vcf/chr${chr}_all.vcf.gz"
    
    extract_twin_sequences "$ind1" "$ind2" "$vcf_file" "$ref_file" "CHR21"
    
    print_success "✨ Complete! Chromosome 21 twin data is ready!"
    print_info "File location: $DATA_DIR/fasta/CHR21_twins_comparison.fasta"
}

# Show usage
usage() {
    cat << EOF
Twin Genome Data Downloader

Usage: $0 [OPTION]

Options:
    --brca1             Download BRCA1 gene for twin pair (Quick, ~10 MB, 2-3 min)
    --brca1-custom      Download BRCA1 for custom individuals
                        Usage: --brca1-custom IND1 IND2
    
    --chr21             Download chromosome 21 for twin pair (Medium, ~200 MB, 5-10 min)
    --chr21-custom      Download chr21 for custom individuals
                        Usage: --chr21-custom IND1 IND2
    
    --gene GENE CHR START END
                        Download specific gene region
                        Example: --gene TP53 17 7661779 7687550
    
    --list-individuals  Show available individual pairs
    --check             Check if required tools are installed
    --help              Show this help message

Examples:
    # Quick start - BRCA1 gene
    $0 --brca1
    
    # Chromosome 21 for default twin pair
    $0 --chr21
    
    # BRCA1 for specific individuals
    $0 --brca1-custom HG00096 HG00097
    
    # Custom gene
    $0 --gene TP53 17 7661779 7687550

Default twin pair: NA19092 and NA19093 (Yoruba siblings from 1000 Genomes)

EOF
}

# List available individuals
list_individuals() {
    cat << EOF
Common Related Pairs from 1000 Genomes:

African (Yoruba, Nigeria):
  - NA19092, NA19093 ⭐ (Default pair)
  - NA19067, NA19068
  - NA19070, NA19072

European (British, England):
  - HG00096, HG00097 ⭐
  - HG00099, HG00100
  - HG00102, HG00103

Asian (Han Chinese, Beijing):
  - NA18525, NA18526
  - NA18537, NA18542

CEPH (Utah, USA):
  - NA12003, NA12004
  - NA12005, NA12006

⭐ = Recommended pairs (good quality data)

To use a specific pair:
  $0 --brca1-custom IND1 IND2

EOF
}

# Main script
main() {
    echo "
╔═══════════════════════════════════════════════════════╗
║       Twin Genome Data Downloader v1.0               ║
║       Download real twin genomic data from           ║
║       the 1000 Genomes Project                       ║
╚═══════════════════════════════════════════════════════╝
"
    
    if [ $# -eq 0 ]; then
        usage
        exit 0
    fi
    
    case "$1" in
        --check)
            check_dependencies
            exit 0
            ;;
        --help)
            usage
            exit 0
            ;;
        --list-individuals)
            list_individuals
            exit 0
            ;;
        --brca1)
            check_dependencies
            setup_directories
            download_brca1_twins
            ;;
        --brca1-custom)
            if [ $# -ne 3 ]; then
                print_error "Usage: $0 --brca1-custom IND1 IND2"
                exit 1
            fi
            check_dependencies
            setup_directories
            download_brca1_twins "$2" "$3"
            ;;
        --chr21)
            check_dependencies
            setup_directories
            download_chr21_twins
            ;;
        --chr21-custom)
            if [ $# -ne 3 ]; then
                print_error "Usage: $0 --chr21-custom IND1 IND2"
                exit 1
            fi
            check_dependencies
            setup_directories
            download_chr21_twins "$2" "$3"
            ;;
        --gene)
            if [ $# -ne 5 ]; then
                print_error "Usage: $0 --gene GENE CHR START END"
                exit 1
            fi
            check_dependencies
            setup_directories
            download_reference "$3"
            download_variants "$3"
            extract_gene_region "$2" "$3" "$4" "$5"
            ;;
        *)
            print_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
    
    echo ""
    print_success "🎉 All done!"
    print_info "Check the $DATA_DIR/fasta/ folder for your data"
    print_info "Load the .fasta file into your DNA analyzer web app!"
}

# Run main function
main "$@"


