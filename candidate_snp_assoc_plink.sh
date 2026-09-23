#!/bin/bash

set -uo pipefail

#===============================================================================
# HELP FUNCTION
#===============================================================================
show_help() {
    cat << HELPEOF

================================================================================
CANDIDATE SNP ASSOCIATION ANALYSIS PIPELINE
================================================================================

DESCRIPTION:
    This script performs genetic association testing for user-specified candidate
    SNPs against a phenotype of interest. It handles PCA calculation for population
    stratification, supports multiple genetic models, and produces formatted
    summary tables.


OPTIONS:
    -h, --help              Show this help message and exit
    -s, --snp-file FILE     Path to SNP definition file (required)
    -g, --geno PREFIX       Path/prefix to PLINK bfiles (required)
    -p, --pheno FILE        Path to phenotype file (required)
    -c, --covar FILE        Path to covariate file (required)
    -o, --output NAME       Output name prefix (required)
    -d, --outdir DIR        Output directory (default: current directory)
    -b, --build BUILD       Genome build: hg19 or hg38 (default: hg19)
    -m, --model MODEL       Genetic model(s): additive, recessive, dominant
                            Can specify multiple separated by comma (default: additive)
    --pheno-name NAME       Phenotype name for labeling (default: PHENO)
    --pheno-col COL         Column name in phenotype file (default: PHENO)
    --covar-names NAMES     Comma-separated covariate names (default: AGE,SEX)
    --n-pcs N               Number of PCs to calculate (default: 10)
    --ld-window N           LD pruning window size (default: 50)
    --ld-step N             LD pruning step size (default: 10)
    --ld-r2 R2              LD pruning r2 threshold (default: 0.2)

SNP FILE FORMAT:
    Tab-delimited file with header. Required columns:
        SNP_ID      - SNP identifier (rsID or custom name)
        CHR         - Chromosome number (1-22, X, Y)
        POS_HG19    - Position in hg19 coordinates
        POS_HG38    - Position in hg38 coordinates
    
    Optional columns:
        GENE        - Gene name (for annotation)
        GROUP       - Group label (for organizing results)

    Example:
        SNP_ID          CHR     POS_HG19    POS_HG38    GENE    GROUP
        rs73885319      22      36661906    36265860    APOL1   G1
        rs60910145      22      36662034    36265988    APOL1   G1
        rs71785313      22      36662046    36266000    APOL1   G2
        rs11912763      22      36684722    36288676    MYH9    MYH9
        rs5750248       22      36702892    36306846    MYH9    MYH9

GENETIC MODELS:
    additive    - Tests additive effect (0, 1, 2 copies of effect allele)
    recessive   - Tests recessive effect (0+1 vs 2 copies)
    dominant    - Tests dominant effect (0 vs 1+2 copies)
    
    Multiple models can be specified: --model additive,recessive,dominant

OUTPUT FILES:
    <output>_summary.txt              - Main summary file with all results
    <output>_results_table.txt        - Formatted results table
    <output>_covariates_with_PCs.txt  - Covariate file with calculated PCs
    <output>_PCs.txt                  - Principal components
    <output>_eigenvalues.txt          - PCA eigenvalues and variance explained
    <output>_<model>.*.glm.linear     - Raw PLINK2 output files
    <output>_<model>.afreq            - Allele frequency files

EXAMPLES:
    # Basic usage with defaults
    sbatch $(basename $0) -s snps.txt -g geno_prefix -p pheno.txt -c covar.txt -o my_analysis

    # Specify multiple models and hg38 build
    sbatch $(basename $0) -s snps.txt -g geno_prefix -p pheno.txt -c covar.txt \\
        -o my_analysis -b hg38 -m additive,recessive,dominant

    # Full specification
    sbatch $(basename $0) -s snps.txt -g geno_prefix -p pheno.txt -c covar.txt \\
        -o kidney_gwas -d /path/to/output -b hg19 -m recessive \\
        --pheno-name eGFR --pheno-col eGFR_normalized --covar-names AGE,SEX,BMI \\
        --n-pcs 10

NOTES:
    - SNPs not found in the genotype data will be skipped with a warning
    - The script automatically calculates PCs for population stratification
    - High-LD regions (MHC, etc.) are excluded from PCA calculation
    - Bonferroni correction is calculated based on total number of tests

DEPENDENCIES:
    - PLINK2 (plink/2.00a2.3 or compatible)
    - R with data.table package

AUTHOR:
    Genetic Analysis Pipeline

================================================================================
HELPEOF
    exit 0
}

#===============================================================================
# PARSE COMMAND LINE ARGUMENTS
#===============================================================================
# Default values
SNP_FILE=""
GENO_PREFIX=""
PHENO_FILE=""
COVAR_FILE=""
OUTPUT_NAME=""
OUTPUT_DIR="."
GENOME_BUILD="hg19"
MODELS="additive"
PHENO_NAME="PHENO"
PHENO_COL="PHENO"
BASE_COVARS="AGE,SEX"
N_PCS=10
LD_WINDOW=50
LD_STEP=10
LD_R2=0.2
USE_CLI=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            ;;
        -s|--snp-file)
            SNP_FILE="$2"
            USE_CLI=true
            shift 2
            ;;
        -g|--geno)
            GENO_PREFIX="$2"
            USE_CLI=true
            shift 2
            ;;
        -p|--pheno)
            PHENO_FILE="$2"
            USE_CLI=true
            shift 2
            ;;
        -c|--covar)
            COVAR_FILE="$2"
            USE_CLI=true
            shift 2
            ;;
        -o|--output)
            OUTPUT_NAME="$2"
            USE_CLI=true
            shift 2
            ;;
        -d|--outdir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -b|--build)
            GENOME_BUILD="$2"
            shift 2
            ;;
        -m|--model)
            MODELS="$2"
            shift 2
            ;;
        --pheno-name)
            PHENO_NAME="$2"
            shift 2
            ;;
        --pheno-col)
            PHENO_COL="$2"
            shift 2
            ;;
        --covar-names)
            BASE_COVARS="$2"
            shift 2
            ;;
        --n-pcs)
            N_PCS="$2"
            shift 2
            ;;
        --ld-window)
            LD_WINDOW="$2"
            shift 2
            ;;
        --ld-step)
            LD_STEP="$2"
            shift 2
            ;;
        --ld-r2)
            LD_R2="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

#===============================================================================
# USER-DEFINED INPUTS (used if not running with command line arguments)
#===============================================================================
if [ "$USE_CLI" = false ]; then
    
    # Path to SNP definition file (see --help for format)
    SNP_FILE="/path/to/snp_list.txt"
    
    # Genotype file prefix (PLINK bfile format)
    GENO_PREFIX="/path/to/genotype_prefix"
    
    # Phenotype file
    PHENO_FILE="/path/to/phenotype.txt"
    
    # Covariate file
    COVAR_FILE="/path/to/covariates.txt"
    
    # Output name prefix
    OUTPUT_NAME="my_analysis"
    
    # Output directory
    OUTPUT_DIR="/path/to/output"
    
    # Genome build: hg19 or hg38
    GENOME_BUILD="hg19"
    
    # Genetic models to test: additive, recessive, dominant (comma-separated)
    MODELS="additive,recessive"
    
    # Phenotype settings
    PHENO_NAME="uPCR"
    PHENO_COL="uPCR_normalized"
    
    # Covariates
    BASE_COVARS="AGE,SEX"
    
    # PCA settings
    N_PCS=10
    LD_WINDOW=50
    LD_STEP=10
    LD_R2=0.2
fi

#===============================================================================
# VALIDATE INPUTS
#===============================================================================
echo ""
echo "========================================================================"
echo "CANDIDATE SNP ASSOCIATION ANALYSIS"
echo "========================================================================"
echo "Date: $(date)"
echo ""

# Check required inputs
missing_inputs=false

if [ -z "$SNP_FILE" ] || [ ! -f "$SNP_FILE" ]; then
    echo "ERROR: SNP file not found: $SNP_FILE"
    missing_inputs=true
fi

if [ -z "$GENO_PREFIX" ] || [ ! -f "${GENO_PREFIX}.bim" ]; then
    echo "ERROR: Genotype files not found: ${GENO_PREFIX}.bim/bed/fam"
    missing_inputs=true
fi

if [ -z "$PHENO_FILE" ] || [ ! -f "$PHENO_FILE" ]; then
    echo "ERROR: Phenotype file not found: $PHENO_FILE"
    missing_inputs=true
fi

if [ -z "$COVAR_FILE" ] || [ ! -f "$COVAR_FILE" ]; then
    echo "ERROR: Covariate file not found: $COVAR_FILE"
    missing_inputs=true
fi

if [ -z "$OUTPUT_NAME" ]; then
    echo "ERROR: Output name not specified"
    missing_inputs=true
fi

if [ "$missing_inputs" = true ]; then
    echo ""
    echo "Use --help for usage information"
    exit 1
fi

if [ "${GENOME_BUILD}" != "hg19" ] && [ "${GENOME_BUILD}" != "hg38" ]; then
    echo "ERROR: GENOME_BUILD must be hg19 or hg38 (got: ${GENOME_BUILD})"
    exit 1
fi

# Validate models
IFS=',' read -ra MODEL_ARRAY <<< "$MODELS"
for model in "${MODEL_ARRAY[@]}"; do
    if [ "$model" != "additive" ] && [ "$model" != "recessive" ] && [ "$model" != "dominant" ]; then
        echo "ERROR: Invalid model '$model'. Must be: additive, recessive, or dominant"
        exit 1
    fi
done

#===============================================================================
# SETUP
#===============================================================================
module load plink/2.00a2.3
module load r

RESULTS_DIR="${OUTPUT_DIR}/${OUTPUT_NAME}"
SNPLIST_DIR="${RESULTS_DIR}/snplists"
TEMP_DIR="${RESULTS_DIR}/temp_pca_$$"

mkdir -p ${RESULTS_DIR} ${SNPLIST_DIR} ${TEMP_DIR}

echo "Configuration:"
echo "  SNP file:      $SNP_FILE"
echo "  Genotype:      $GENO_PREFIX"
echo "  Phenotype:     $PHENO_NAME ($PHENO_COL)"
echo "  Covariates:    $BASE_COVARS + ${N_PCS} PCs"
echo "  Genome build:  $GENOME_BUILD"
echo "  Models:        $MODELS"
echo "  Output:        $RESULTS_DIR"
echo ""

# Initialize summary file
SUMMARY_FILE="${RESULTS_DIR}/${OUTPUT_NAME}_summary.txt"
cat > ${SUMMARY_FILE} << EOF
================================================================================
CANDIDATE SNP ASSOCIATION ANALYSIS SUMMARY
================================================================================
Date: $(date)
Output name: ${OUTPUT_NAME}

INPUT FILES
-----------
SNP file: ${SNP_FILE}
Genotype: ${GENO_PREFIX}
Phenotype file: ${PHENO_FILE}
Covariate file: ${COVAR_FILE}

PARAMETERS
----------
Phenotype: ${PHENO_NAME} (column: ${PHENO_COL})
Genome build: ${GENOME_BUILD}
Genetic models: ${MODELS}
Base covariates: ${BASE_COVARS}
Number of PCs: ${N_PCS}
LD pruning: window=${LD_WINDOW}, step=${LD_STEP}, r2=${LD_R2}

EOF

#===============================================================================
# STEP 1: PCA CALCULATION
#===============================================================================
echo "========================================================================"
echo "STEP 1: Principal Component Analysis"
echo "========================================================================"

cat > ${TEMP_DIR}/high_ld_regions.txt << EOF
5 44000000 51500000 r1
6 25000000 33500000 r2
8 8000000 12000000 r3
11 45000000 57000000 r4
EOF

echo "LD pruning..."
plink2 --bfile ${GENO_PREFIX} \
    --exclude range ${TEMP_DIR}/high_ld_regions.txt \
    --indep-pairwise ${LD_WINDOW} ${LD_STEP} ${LD_R2} \
    --out ${TEMP_DIR}/pruned_snps \
    --silent

if [ ! -f ${TEMP_DIR}/pruned_snps.prune.in ]; then
    echo "ERROR: LD pruning failed"
    exit 1
fi

N_PRUNED=$(wc -l < ${TEMP_DIR}/pruned_snps.prune.in)
echo "  SNPs after pruning: ${N_PRUNED}"

echo "Calculating PCs..."
plink2 --bfile ${GENO_PREFIX} \
    --extract ${TEMP_DIR}/pruned_snps.prune.in \
    --pca ${N_PCS} \
    --out ${TEMP_DIR}/pca \
    --silent

if [ ! -f ${TEMP_DIR}/pca.eigenvec ]; then
    echo "ERROR: PCA calculation failed"
    exit 1
fi

echo "Merging covariates with PCs..."
Rscript --vanilla - "${TEMP_DIR}" "${COVAR_FILE}" "${N_PCS}" "${RESULTS_DIR}" "${OUTPUT_NAME}" "${BASE_COVARS}" << 'REOF'
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
temp_dir <- args[1]; covar_file <- args[2]; n_pcs <- as.integer(args[3])
results_dir <- args[4]; output_name <- args[5]; base_covars <- args[6]

pca <- fread(file.path(temp_dir, "pca.eigenvec"))
if(names(pca)[1]=="#FID") names(pca)[1] <- "FID"
pc_names <- paste0("PC", 1:n_pcs)
names(pca)[3:(2+n_pcs)] <- pc_names

covar <- fread(covar_file, na.strings=c("","NA","N/A","."))
if(!("FID" %in% names(covar))) {
    if("IID" %in% names(covar)) { 
        names(covar)[1] <- "FID"
        if(names(covar)[2]!="IID") names(covar)[2] <- "IID" 
    } else {
        names(covar)[1:2] <- c("FID","IID")
    }
}

base_list <- unlist(strsplit(base_covars, ","))
available <- intersect(base_list, names(covar))
keep_cols <- c("FID", "IID", available)
covar_subset <- covar[, ..keep_cols]

merged <- merge(covar_subset, pca[, c("FID","IID",pc_names), with=FALSE], by=c("FID","IID"))
all_covar_cols <- c(available, pc_names)
merged_complete <- merged[complete.cases(merged[, ..all_covar_cols])]

final_cols <- c("FID", "IID", available, pc_names)
merged_complete <- merged_complete[, ..final_cols]

output_file <- file.path(results_dir, paste0(output_name, "_covariates_with_PCs.txt"))
write.table(merged_complete, output_file, sep="\t", row.names=FALSE, quote=FALSE, na="")

fwrite(pca[, c("FID","IID",pc_names), with=FALSE], 
       file.path(results_dir, paste0(output_name,"_PCs.txt")), sep="\t")

if(file.exists(file.path(temp_dir,"pca.eigenval"))) {
    ev <- fread(file.path(temp_dir,"pca.eigenval"), header=FALSE)
    ev$PC <- paste0("PC", 1:nrow(ev))
    ev$pct <- round(ev$V1/sum(ev$V1)*100, 2)
    fwrite(ev[, .(PC, eigenvalue=V1, variance_pct=pct)], 
           file.path(results_dir, paste0(output_name,"_eigenvalues.txt")), sep="\t")
    cat("Variance explained by top 5 PCs:", sum(ev$pct[1:min(5,nrow(ev))]), "%\n")
}

cat("Samples with complete data:", nrow(merged_complete), "\n")
REOF

if [ ! -f ${RESULTS_DIR}/${OUTPUT_NAME}_covariates_with_PCs.txt ]; then
    echo "ERROR: Failed to create covariate file"
    exit 1
fi

COVAR_FILE_WITH_PCS="${RESULTS_DIR}/${OUTPUT_NAME}_covariates_with_PCs.txt"
PC_NAMES=$(seq -s "," -f "PC%.0f" 1 ${N_PCS})
COVAR_NAMES="${BASE_COVARS},${PC_NAMES}"

echo "  PCA complete"
echo ""

cat >> ${SUMMARY_FILE} << EOF

PCA RESULTS
-----------
SNPs used for PCA: ${N_PRUNED}
PCs calculated: ${N_PCS}

EOF

#===============================================================================
# STEP 2: PARSE SNP FILE AND IDENTIFY SNPS
#===============================================================================
echo "========================================================================"
echo "STEP 2: SNP Identification"
echo "========================================================================"

BIM_FILE="${GENO_PREFIX}.bim"

# Use R to parse SNP file and find matches
Rscript --vanilla - "${SNP_FILE}" "${BIM_FILE}" "${GENOME_BUILD}" "${SNPLIST_DIR}" "${RESULTS_DIR}" "${OUTPUT_NAME}" << 'REOF'
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
snp_file <- args[1]
bim_file <- args[2]
genome_build <- args[3]
snplist_dir <- args[4]
results_dir <- args[5]
output_name <- args[6]

cat("Reading SNP definition file...\n")
snps <- fread(snp_file)

# Standardize column names
names(snps) <- toupper(names(snps))
if("RSID" %in% names(snps) && !("SNP_ID" %in% names(snps))) {
    names(snps)[names(snps)=="RSID"] <- "SNP_ID"
}
if("ID" %in% names(snps) && !("SNP_ID" %in% names(snps))) {
    names(snps)[names(snps)=="ID"] <- "SNP_ID"
}

required_cols <- c("SNP_ID", "CHR")
pos_col <- ifelse(genome_build == "hg19", "POS_HG19", "POS_HG38")

if(!all(required_cols %in% names(snps))) {
    stop("SNP file must contain columns: SNP_ID, CHR")
}
if(!(pos_col %in% names(snps))) {
    stop(paste("SNP file must contain column:", pos_col))
}

cat("SNPs in input file:", nrow(snps), "\n")
cat("Using position column:", pos_col, "\n")

# Read BIM file
cat("Reading BIM file...\n")
bim <- fread(bim_file, header=FALSE, col.names=c("CHR","ID","CM","POS","A1","A2"))
cat("Variants in genotype data:", nrow(bim), "\n")

# Match SNPs by chromosome and position
snps$POS <- snps[[pos_col]]
snps$CHR <- as.character(snps$CHR)
bim$CHR <- as.character(bim$CHR)

snps$FOUND_ID <- NA_character_
snps$FOUND <- FALSE

for(i in 1:nrow(snps)) {
    chr <- snps$CHR[i]
    pos <- snps$POS[i]
    
    match_idx <- which(bim$CHR == chr & bim$POS == pos)
    
    if(length(match_idx) > 0) {
        snps$FOUND_ID[i] <- bim$ID[match_idx[1]]
        snps$FOUND[i] <- TRUE
    }
}

cat("\nSNP matching results:\n")
cat("  Found:", sum(snps$FOUND), "\n")
cat("  Not found:", sum(!snps$FOUND), "\n\n")

# Report found/not found
if(any(snps$FOUND)) {
    cat("Found SNPs:\n")
    found_snps <- snps[FOUND == TRUE]
    for(i in 1:nrow(found_snps)) {
        cat("  ", found_snps$SNP_ID[i], " (chr", found_snps$CHR[i], ":", 
            found_snps$POS[i], ") -> ", found_snps$FOUND_ID[i], "\n", sep="")
    }
}

if(any(!snps$FOUND)) {
    cat("\nNot found SNPs:\n")
    notfound_snps <- snps[FOUND == FALSE]
    for(i in 1:nrow(notfound_snps)) {
        cat("  ", notfound_snps$SNP_ID[i], " (chr", notfound_snps$CHR[i], ":", 
            notfound_snps$POS[i], ")\n", sep="")
    }
}

# Write SNP list for PLINK
found_snps <- snps[FOUND == TRUE]
if(nrow(found_snps) == 0) {
    stop("ERROR: No SNPs found in genotype data!")
}

writeLines(found_snps$FOUND_ID, file.path(snplist_dir, "all_snps.txt"))

# Save SNP mapping info
snps$GENOME_BUILD <- genome_build
fwrite(snps, file.path(results_dir, paste0(output_name, "_snp_mapping.txt")), sep="\t")

# Save found SNPs info for later use
fwrite(found_snps[, .(SNP_ID, CHR, POS, FOUND_ID, 
                       GENE = ifelse("GENE" %in% names(found_snps), get("GENE"), NA),
                       GROUP = ifelse("GROUP" %in% names(found_snps), get("GROUP"), NA))],
       file.path(snplist_dir, "found_snps_info.txt"), sep="\t")

cat("\nSNP list written to:", file.path(snplist_dir, "all_snps.txt"), "\n")
REOF

if [ ! -f ${SNPLIST_DIR}/all_snps.txt ]; then
    echo "ERROR: SNP identification failed"
    exit 1
fi

N_SNPS=$(wc -l < ${SNPLIST_DIR}/all_snps.txt)
echo "  SNPs to test: ${N_SNPS}"
echo ""

cat >> ${SUMMARY_FILE} << EOF

SNP IDENTIFICATION
------------------
EOF
cat ${RESULTS_DIR}/${OUTPUT_NAME}_snp_mapping.txt >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}

#===============================================================================
# STEP 3: ASSOCIATION TESTING
#===============================================================================
echo "========================================================================"
echo "STEP 3: Association Testing"
echo "========================================================================"

TOTAL_TESTS=0

for model in "${MODEL_ARRAY[@]}"; do
    echo ""
    echo "Running ${model} model..."
    
    case $model in
        additive)
            plink2 --bfile ${GENO_PREFIX} \
                --extract ${SNPLIST_DIR}/all_snps.txt \
                --pheno ${PHENO_FILE} --pheno-name ${PHENO_COL} \
                --covar ${COVAR_FILE_WITH_PCS} --covar-name ${COVAR_NAMES} \
                --glm hide-covar --freq \
                --out ${RESULTS_DIR}/${OUTPUT_NAME}_additive
            ;;
        recessive)
            plink2 --bfile ${GENO_PREFIX} \
                --extract ${SNPLIST_DIR}/all_snps.txt \
                --pheno ${PHENO_FILE} --pheno-name ${PHENO_COL} \
                --covar ${COVAR_FILE_WITH_PCS} --covar-name ${COVAR_NAMES} \
                --glm recessive hide-covar --freq \
                --out ${RESULTS_DIR}/${OUTPUT_NAME}_recessive
            ;;
        dominant)
            plink2 --bfile ${GENO_PREFIX} \
                --extract ${SNPLIST_DIR}/all_snps.txt \
                --pheno ${PHENO_FILE} --pheno-name ${PHENO_COL} \
                --covar ${COVAR_FILE_WITH_PCS} --covar-name ${COVAR_NAMES} \
                --glm dominant hide-covar --freq \
                --out ${RESULTS_DIR}/${OUTPUT_NAME}_dominant
            ;;
    esac
    
    # Check for output
    result_file="${RESULTS_DIR}/${OUTPUT_NAME}_${model}.${PHENO_COL}.glm.linear"
    if [ -f "$result_file" ]; then
        n_results=$(tail -n +2 "$result_file" | wc -l)
        echo "  ${model}: ${n_results} tests completed"
        TOTAL_TESTS=$((TOTAL_TESTS + n_results))
    else
        echo "  ${model}: No results (check for errors)"
    fi
done

echo ""
echo "Total tests: ${TOTAL_TESTS}"

#===============================================================================
# STEP 4: CREATE SUMMARY TABLE
#===============================================================================
echo ""
echo "========================================================================"
echo "STEP 4: Creating Summary Table"
echo "========================================================================"

Rscript --vanilla - "${RESULTS_DIR}" "${OUTPUT_NAME}" "${MODELS}" "${PHENO_COL}" "${SNPLIST_DIR}" "${N_SNPS}" << 'REOF'
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
results_dir <- args[1]
output_name <- args[2]
models <- unlist(strsplit(args[3], ","))
pheno_col <- args[4]
snplist_dir <- args[5]
n_snps <- as.integer(args[6])

cat("Creating summary table...\n")

# Read SNP info
snp_info <- fread(file.path(snplist_dir, "found_snps_info.txt"))

# Collect all results
all_results <- list()

for(model in models) {
    result_file <- file.path(results_dir, paste0(output_name, "_", model, ".", pheno_col, ".glm.linear"))
    freq_file <- file.path(results_dir, paste0(output_name, "_", model, ".afreq"))
    
    if(file.exists(result_file)) {
        res <- fread(result_file)
        if(names(res)[1] == "#CHROM") names(res)[1] <- "CHROM"
        
        # Filter to actual test results (not NA rows)
        res <- res[!is.na(P)]
        
        if(nrow(res) > 0) {
            res$MODEL <- toupper(model)
            
            # Add frequency info if available
            if(file.exists(freq_file)) {
                freq <- fread(freq_file)
                if(names(freq)[1] == "#CHROM") names(freq)[1] <- "CHROM"
                res <- merge(res, freq[, .(ID, ALT_FREQS)], by="ID", all.x=TRUE)
                names(res)[names(res)=="ALT_FREQS"] <- "MAF"
            }
            
            all_results[[model]] <- res
        }
    }
}

if(length(all_results) == 0) {
    cat("WARNING: No results to summarize\n")
    quit(status=0)
}

# Combine results
combined <- rbindlist(all_results, fill=TRUE)

# Merge with SNP info
combined <- merge(combined, snp_info[, .(FOUND_ID, SNP_ID, GENE, GROUP)], 
                  by.x="ID", by.y="FOUND_ID", all.x=TRUE)

# Calculate Bonferroni threshold
n_tests <- nrow(combined)
bonf_thresh <- 0.05 / n_tests

# Add significance flags
combined$BONF_SIG <- combined$P < bonf_thresh
combined$NOM_SIG <- combined$P < 0.05

# Select and order columns for output
out_cols <- c("SNP_ID", "ID", "CHROM", "POS", "GENE", "GROUP", "MODEL",
              "A1", "REF", "ALT", "MAF", "OBS_CT", "BETA", "SE", "T_STAT", "P",
              "NOM_SIG", "BONF_SIG")
out_cols <- intersect(out_cols, names(combined))
combined <- combined[, ..out_cols]

# Sort by P-value
combined <- combined[order(P)]

# Save full results table
fwrite(combined, file.path(results_dir, paste0(output_name, "_results_table.txt")), sep="\t")

# Create formatted summary
cat("\n")
cat("================================================================================\n")
cat("ASSOCIATION RESULTS SUMMARY\n")
cat("================================================================================\n")
cat("\n")
cat("Total tests:", n_tests, "\n")
cat("Bonferroni threshold:", format(bonf_thresh, scientific=TRUE, digits=3), "\n")
cat("Nominally significant (P<0.05):", sum(combined$NOM_SIG, na.rm=TRUE), "\n")
cat("Bonferroni significant:", sum(combined$BONF_SIG, na.rm=TRUE), "\n")
cat("\n")

# Print top results
cat("TOP RESULTS (sorted by P-value):\n")
cat(paste(rep("-", 80), collapse=""), "\n")

print_cols <- c("SNP_ID", "GENE", "MODEL", "BETA", "SE", "P", "MAF")
print_cols <- intersect(print_cols, names(combined))
top_n <- min(20, nrow(combined))

print(combined[1:top_n, ..print_cols], digits=3)

# Results by model
cat("\n")
cat("RESULTS BY MODEL:\n")
cat(paste(rep("-", 80), collapse=""), "\n")
for(m in models) {
    m_upper <- toupper(m)
    m_res <- combined[MODEL == m_upper]
    if(nrow(m_res) > 0) {
        cat("\n", m_upper, " (n=", nrow(m_res), "):\n", sep="")
        cat("  Min P-value:", format(min(m_res$P, na.rm=TRUE), scientific=TRUE, digits=3), "\n")
        cat("  Nominally significant:", sum(m_res$NOM_SIG, na.rm=TRUE), "\n")
        cat("  Bonferroni significant:", sum(m_res$BONF_SIG, na.rm=TRUE), "\n")
        
        # Top hit for this model
        top_hit <- m_res[which.min(P)]
        if(nrow(top_hit) > 0) {
            cat("  Top hit:", top_hit$SNP_ID, 
                "(BETA=", round(top_hit$BETA, 4), 
                ", P=", format(top_hit$P, scientific=TRUE, digits=3), ")\n")
        }
    }
}

# Save summary stats
summary_stats <- data.table(
    Metric = c("Total_tests", "Bonferroni_threshold", "Nominally_significant", "Bonferroni_significant"),
    Value = c(n_tests, bonf_thresh, sum(combined$NOM_SIG, na.rm=TRUE), sum(combined$BONF_SIG, na.rm=TRUE))
)
fwrite(summary_stats, file.path(results_dir, paste0(output_name, "_summary_stats.txt")), sep="\t")

cat("\n")
cat("Results table saved to:", file.path(results_dir, paste0(output_name, "_results_table.txt")), "\n")
REOF

#===============================================================================
# STEP 5: FINALIZE
#===============================================================================
echo ""
echo "========================================================================"
echo "FINALIZING"
echo "========================================================================"

# Calculate Bonferroni threshold
if [ "$TOTAL_TESTS" -gt 0 ]; then
    BONF=$(echo "scale=6; 0.05 / ${TOTAL_TESTS}" | bc)
else
    BONF="NA"
fi

# Append final summary
cat >> ${SUMMARY_FILE} << EOF

================================================================================
RESULTS SUMMARY
================================================================================
Total association tests: ${TOTAL_TESTS}
Bonferroni threshold (alpha=0.05): ${BONF}

Models tested: ${MODELS}
SNPs tested: ${N_SNPS}

OUTPUT FILES
------------
EOF

ls -1 ${RESULTS_DIR}/${OUTPUT_NAME}* >> ${SUMMARY_FILE} 2>/dev/null

cat >> ${SUMMARY_FILE} << EOF

Analysis completed: $(date)
================================================================================
EOF

# Cleanup
rm -rf ${TEMP_DIR}

echo ""
echo "Analysis complete!"
echo ""
echo "Output directory: ${RESULTS_DIR}"
echo "Summary file: ${SUMMARY_FILE}"
echo "Results table: ${RESULTS_DIR}/${OUTPUT_NAME}_results_table.txt"
echo ""
echo "Key output files:"
ls -lh ${RESULTS_DIR}/${OUTPUT_NAME}_results_table.txt 2>/dev/null
ls -lh ${RESULTS_DIR}/${OUTPUT_NAME}_summary_stats.txt 2>/dev/null
echo ""
echo "Tests: ${TOTAL_TESTS}, Bonferroni threshold: ${BONF}"
echo ""
echo "Done! $(date)"
