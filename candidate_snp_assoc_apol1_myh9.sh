#!/bin/bash

#SBATCH --account [USER_DEFINED]
#SBATCH --partition [USER_DEFINED]
#SBATCH --time [USER_DEFINED]
#SBATCH --mem [USER_DEFINED]
#SBATCH --nodes [USER_DEFINED]
#SBATCH --ntasks [USER_DEFINED]
#SBATCH --job-name candidate-snp-assoc
#SBATCH --output slurm_%j.log

set -uo pipefail

module load plink/2.00a2.3
module load r

#################### USER-DEFINED INPUTS ###########################

COHORT="[USER_DEFINED]" # cohort identifier (will be used in output file naming and in summary table creation)
PHENO_NAME="[USER_DEFINED]" # phenotype name (will be used in output file naming and in summary table creation)
PHENO_COL="[USER_DEFINED]" #name of column containing trait for association in phenotype file
GENOME_BUILD="[USER_DEFINED]" #[hg19/hg38]

GENO_PREFIX="[USER_DEFINED]" #full path to plink file 
PHENO_FILE="[USER_DEFINED]" #full path to phenotype file 
COVAR_FILE="[USER_DEFINED]" #full path to covariate file 

BASE_COVARS="[USER_DEFINED]" # comma-seperated list of covariates as matched with covariate file
N_PCS=[USER_DEFINED] #number of PCs
LD_WINDOW= [USER_DEFINED] #number
LD_STEP=[USER_DEFINED] #number
LD_R2=[USER_DEFINED] #number

OUTPUT_DIR="/data/awonkam1/fujr/SCD_Kidney/candidate_SNP_assoc"

####################################################################

OUTPUT_PREFIX="$${COHORT}_$${PHENO_NAME}"
RESULTS_DIR="$${OUTPUT_DIR}/$${COHORT}/${PHENO_NAME}"
SNPLIST_DIR="${RESULTS_DIR}/snplists"
TEMP_DIR="${RESULTS_DIR}/temp_pca_$$"

mkdir -p $${RESULTS_DIR} $${SNPLIST_DIR} ${TEMP_DIR}

echo "Candidate SNP Association Analysis"
echo "Date: $(date)"
echo "Cohort: ${COHORT}"
echo "Phenotype: $${PHENO_NAME} ($${PHENO_COL})"
echo "Build: ${GENOME_BUILD}"
echo "Covariates: $${BASE_COVARS} + $${N_PCS} PCs"
echo ""

if [ "$${GENOME_BUILD}" != "hg19" ] && [ "$${GENOME_BUILD}" != "hg38" ]; then
    echo "ERROR: GENOME_BUILD must be hg19 or hg38"
    exit 1
fi

SUMMARY_FILE="$${RESULTS_DIR}/$${OUTPUT_PREFIX}_summary.txt"
cat > ${SUMMARY_FILE} << EOF
Candidate SNP Association Results
Date: $(date)
Cohort: ${COHORT}
Phenotype: $${PHENO_NAME} ($${PHENO_COL})
Build: ${GENOME_BUILD}
Genotype: ${GENO_PREFIX}
Covariates: $${BASE_COVARS} + $${N_PCS} PCs
LD pruning: $${LD_WINDOW} $${LD_STEP} ${LD_R2}

EOF

# PCA
echo "Running PCA..."

if [ "${GENOME_BUILD}" == "hg19" ]; then
    cat > ${TEMP_DIR}/high_ld_regions.txt << EOF
5 44000000 51500000 r1
6 25000000 33500000 r2
8 8000000 12000000 r3
11 45000000 57000000 r4
EOF
else
    cat > ${TEMP_DIR}/high_ld_regions.txt << EOF
5 44000000 51500000 r1
6 25000000 33500000 r2
8 8000000 12000000 r3
11 45000000 57000000 r4
EOF
fi

plink2 --bfile ${GENO_PREFIX} \
    --exclude range ${TEMP_DIR}/high_ld_regions.txt \
    --indep-pairwise $${LD_WINDOW} $${LD_STEP} ${LD_R2} \
    --out ${TEMP_DIR}/pruned_snps

if [ ! -f ${TEMP_DIR}/pruned_snps.prune.in ]; then
    echo "ERROR: LD pruning failed"
    exit 1
fi

N_PRUNED=$$(wc -l < $${TEMP_DIR}/pruned_snps.prune.in)
echo "SNPs for PCA: ${N_PRUNED}"

plink2 --bfile ${GENO_PREFIX} \
    --extract ${TEMP_DIR}/pruned_snps.prune.in \
    --pca ${N_PCS} \
    --out ${TEMP_DIR}/pca

if [ ! -f ${TEMP_DIR}/pca.eigenvec ]; then
    echo "ERROR: PCA failed"
    exit 1
fi

Rscript --vanilla - "$${TEMP_DIR}" "$${COVAR_FILE}" "$${N_PCS}" "$${RESULTS_DIR}" "$${OUTPUT_PREFIX}" "$${BASE_COVARS}" << 'EOF'
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
temp_dir <- args[1]; covar_file <- args[2]; n_pcs <- as.integer(args[3])
results_dir <- args[4]; output_prefix <- args[5]; base_covars <- args[6]

pca <- fread(file.path(temp_dir, "pca.eigenvec"))
if(names(pca)[1]=="#FID") names(pca)[1] <- "FID"
pc_names <- paste0("PC", 1:n_pcs)
names(pca)[3:(2+n_pcs)] <- pc_names

covar <- fread(covar_file, na.strings=c("","NA","N/A","."))
if(!("FID" %in% names(covar))) {
    if("IID" %in% names(covar)) { names(covar)[1] <- "FID"; if(names(covar)[2]!="IID") names(covar)[2] <- "IID" }
    else names(covar)[1:2] <- c("FID","IID")
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
output_file <- file.path(results_dir, paste0(output_prefix, "_covariates_with_PCs.txt"))
write.table(merged_complete, output_file, sep="\t", row.names=FALSE, quote=FALSE, na="")

fwrite(pca[, c("FID","IID",pc_names), with=FALSE], file.path(results_dir, paste0(output_prefix,"_PCs.txt")), sep="\t")
if(file.exists(file.path(temp_dir,"pca.eigenval"))) {
    ev <- fread(file.path(temp_dir,"pca.eigenval"), header=FALSE)
    ev$$PC <- paste0("PC", 1:nrow(ev)); ev$$pct <- round(ev$$V1/sum(ev$$V1)*100, 2)
    fwrite(ev[, .(PC, eigenvalue=V1, variance_pct=pct)], file.path(results_dir, paste0(output_prefix,"_eigenvalues.txt")), sep="\t")
}
EOF

if [ ! -f $${RESULTS_DIR}/$${OUTPUT_PREFIX}_covariates_with_PCs.txt ]; then
    echo "ERROR: Covariate file creation failed"
    exit 1
fi

COVAR_FILE_WITH_PCS="$${RESULTS_DIR}/$${OUTPUT_PREFIX}_covariates_with_PCs.txt"
PC_NAMES=$$(seq -s "," -f "PC%.0f" 1 $${N_PCS})
COVAR_NAMES="$${BASE_COVARS},$${PC_NAMES}"

echo "PCA complete"
echo "" >> ${SUMMARY_FILE}
echo "PCA: $${N_PRUNED} SNPs, $${N_PCS} PCs" >> ${SUMMARY_FILE}

# SNP identification
get_snp_id() {
    awk -v c=\$2 -v p=\$3 '\$1==c && \$4==p {print \$2}' \$1
}

BIM_FILE="${GENO_PREFIX}.bim"

declare -A APOL1_POS MYH9_POS

if [ "${GENOME_BUILD}" == "hg19" ]; then
    APOL1_POS=([G1]=36661906 [G1M]=36662034 [G2]=36662046)
    MYH9_POS=([rs11912763]=36684722 [rs16996648]=36692752 [rs5750248]=36702892 [rs1557529]=36705529 [rs8141189]=36714710 [rs1005570]=36715274 [rs16996672]=36725970)
else
    APOL1_POS=([G1]=36265860 [G1M]=36265988 [G2]=36266000)
    MYH9_POS=([rs11912763]=36288676 [rs16996648]=36296706 [rs5750248]=36306846 [rs1557529]=36309484 [rs8141189]=36318665 [rs1005570]=36319229 [rs16996672]=36329925)
fi

echo "Identifying SNPs (${GENOME_BUILD})..."
> $${SNPLIST_DIR}/APOL1_g1.txt; > $${SNPLIST_DIR}/APOL1_g1m.txt
> $${SNPLIST_DIR}/APOL1_g2.txt; > $${SNPLIST_DIR}/APOL1_all.txt

APOL1_G1_FOUND=false; APOL1_G1M_FOUND=false; APOL1_G2_FOUND=false

for snp in G1 G1M G2; do
    pos=$${APOL1_POS[$$snp]}
    snp_id=$$(get_snp_id $${BIM_FILE} 22 ${pos})
    if [ -n "$snp_id" ]; then
        echo "  APOL1 $${snp}: $${snp_id}"
        case $snp in
            G1) echo "$$snp_id" > $${SNPLIST_DIR}/APOL1_g1.txt; echo "$$snp_id" >> $${SNPLIST_DIR}/APOL1_all.txt; APOL1_G1_FOUND=true ;;
            G1M) echo "$$snp_id" > $${SNPLIST_DIR}/APOL1_g1m.txt; echo "$$snp_id" >> $${SNPLIST_DIR}/APOL1_all.txt; APOL1_G1M_FOUND=true ;;
            G2) echo "$$snp_id" > $${SNPLIST_DIR}/APOL1_g2.txt; echo "$$snp_id" >> $${SNPLIST_DIR}/APOL1_all.txt; APOL1_G2_FOUND=true ;;
        esac
    else
        echo "  APOL1 ${snp}: not found"
    fi
done

> ${SNPLIST_DIR}/MYH9_snps.txt
MYH9_COUNT=0
for snp in "${!MYH9_POS[@]}"; do
    pos=$${MYH9_POS[$$snp]}
    snp_id=$$(get_snp_id $${BIM_FILE} 22 ${pos})
    if [ -n "$snp_id" ]; then
        echo "  MYH9 $${snp}: $${snp_id}"
        echo "$$snp_id" >> $${SNPLIST_DIR}/MYH9_snps.txt
        MYH9_COUNT=$((MYH9_COUNT + 1))
    fi
done

APOL1_G1_COMPLETE=false
[ "$$APOL1_G1_FOUND" = true ] && [ "$$APOL1_G1M_FOUND" = true ] && APOL1_G1_COMPLETE=true

RUN_APOL1_COMBINED=false
[ "$$APOL1_G1_COMPLETE" = true ] || [ "$$APOL1_G2_FOUND" = true ] && RUN_APOL1_COMBINED=true

echo ""
echo "APOL1: G1=$${APOL1_G1_FOUND}, G1M=$${APOL1_G1M_FOUND}, G2=${APOL1_G2_FOUND}"
echo "G1 haplotype complete: ${APOL1_G1_COMPLETE}"
echo "MYH9: ${MYH9_COUNT}/7 SNPs"

cat >> ${SUMMARY_FILE} << EOF

SNP Availability
APOL1 G1 (rs73885319): ${APOL1_G1_FOUND}
APOL1 G1-M (rs60910145): ${APOL1_G1M_FOUND}
APOL1 G2 (rs71785313): ${APOL1_G2_FOUND}
G1 haplotype complete: ${APOL1_G1_COMPLETE}
MYH9: ${MYH9_COUNT}/7

EOF

TOTAL_TESTS=0

# APOL1 individual SNP tests
echo ""
echo "Running APOL1 tests..."

if [ "$$APOL1_G1_FOUND" = true ] && [ -s $${SNPLIST_DIR}/APOL1_g1.txt ]; then
    plink2 --bfile $${GENO_PREFIX} --extract $${SNPLIST_DIR}/APOL1_g1.txt \
        --pheno $${PHENO_FILE} --pheno-name $${PHENO_COL} \
        --covar $${COVAR_FILE_WITH_PCS} --covar-name $${COVAR_NAMES} \
        --glm recessive hide-covar --freq \
        --out $${RESULTS_DIR}/$${OUTPUT_PREFIX}_APOL1_G1_recessive
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    echo "  G1 recessive: done"
fi

if [ "$$APOL1_G1M_FOUND" = true ] && [ -s $${SNPLIST_DIR}/APOL1_g1m.txt ]; then
    plink2 --bfile $${GENO_PREFIX} --extract $${SNPLIST_DIR}/APOL1_g1m.txt \
        --pheno $${PHENO_FILE} --pheno-name $${PHENO_COL} \
        --covar $${COVAR_FILE_WITH_PCS} --covar-name $${COVAR_NAMES} \
        --glm recessive hide-covar --freq \
        --out $${RESULTS_DIR}/$${OUTPUT_PREFIX}_APOL1_G1M_recessive
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    echo "  G1-M recessive: done"
fi

if [ "$$APOL1_G2_FOUND" = true ] && [ -s $${SNPLIST_DIR}/APOL1_g2.txt ]; then
    plink2 --bfile $${GENO_PREFIX} --extract $${SNPLIST_DIR}/APOL1_g2.txt \
        --pheno $${PHENO_FILE} --pheno-name $${PHENO_COL} \
        --covar $${COVAR_FILE_WITH_PCS} --covar-name $${COVAR_NAMES} \
        --glm recessive hide-covar --freq \
        --out $${RESULTS_DIR}/$${OUTPUT_PREFIX}_APOL1_G2_recessive
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    echo "  G2 recessive: done"
fi

# APOL1 combined model
if [ "$RUN_APOL1_COMBINED" = true ]; then
    APOL1_COUNT=$$(wc -l < $${SNPLIST_DIR}/APOL1_all.txt 2>/dev/null || echo 0)
    if [ "$APOL1_COUNT" -ge 1 ]; then
        plink2 --bfile $${GENO_PREFIX} --extract $${SNPLIST_DIR}/APOL1_all.txt \
            --recode A --freq --out ${RESULTS_DIR}/apol1_geno

        if [ -f ${RESULTS_DIR}/apol1_geno.raw ]; then
            [ "${GENOME_BUILD}" == "hg19" ] && { G1_POS=36661906; G1M_POS=36662034; G2_POS=36662046; } || { G1_POS=36265860; G1M_POS=36265988; G2_POS=36266000; }

            Rscript --vanilla - "$${RESULTS_DIR}" "$${OUTPUT_PREFIX}" "$${PHENO_FILE}" "$${COVAR_FILE_WITH_PCS}" "$${COVAR_NAMES}" "$${APOL1_G1_FOUND}" "$${APOL1_G1M_FOUND}" "$${APOL1_G2_FOUND}" "$${PHENO_COL}" "$${G1_POS}" "$${G1M_POS}" "$${G2_POS}" "${APOL1_G1_COMPLETE}" << 'EOF'
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
results_dir <- args[1]; output_prefix <- args[2]; pheno_file <- args[3]
covar_file <- args[4]; covar_names <- args[5]
g1_found <- args[6]=="true"; g1m_found <- args[7]=="true"; g2_found <- args[8]=="true"
pheno_col <- args[9]; g1_pos <- args[10]; g1m_pos <- args[11]; g2_pos <- args[12]
g1_complete <- args[13]=="true"

geno <- fread(file.path(results_dir, "apol1_geno.raw"))
geno$$G1 <- 0; geno$$G1M <- 0; geno$G2 <- 0

if(g1_found) { g1_col <- grep(g1_pos, names(geno), value=TRUE); if(length(g1_col)>0) geno$G1 <- geno[[g1_col[1]]] }
if(g1m_found) { g1m_col <- grep(g1m_pos, names(geno), value=TRUE); if(length(g1m_col)>0) geno$G1M <- geno[[g1m_col[1]]] }
if(g2_found) { g2_col <- grep(g2_pos, names(geno), value=TRUE); if(length(g2_col)>0) geno$G2 <- geno[[g2_col[1]]] }

geno$$G1[is.na(geno$$G1)] <- 0; geno$$G1M[is.na(geno$$G1M)] <- 0; geno$$G2[is.na(geno$$G2)] <- 0

if(g1_complete) { geno$$G1_haplotype <- pmin(geno$$G1, geno$$G1M) } else { geno$$G1_haplotype <- 0 }
geno$$total_risk <- geno$$G1_haplotype + geno$G2
geno$$high_risk <- as.integer(geno$$total_risk >= 2)

fwrite(geno[, .(FID, IID, G1, G1M, G2, G1_haplotype, total_risk, high_risk)],
       file.path(results_dir, paste0(output_prefix, "_APOL1_genotypes.txt")), sep="\t")

n_high <- sum(geno$high_risk==1, na.rm=TRUE)
if(n_high < 3) {
    results <- data.table(Test="APOL1_combined_recessive", Beta=NA, SE=NA, P_value=NA, N=nrow(geno), N_high_risk=n_high, Note="Insufficient high-risk N")
    fwrite(results, file.path(results_dir, paste0(output_prefix, "_APOL1_combined_results.txt")), sep="\t")
    quit(status=0)
}

pheno <- fread(pheno_file); covar <- fread(covar_file)
if(!("FID" %in% names(pheno))) { if("IID" %in% names(pheno)) names(pheno)[1] <- "FID" else names(pheno)[1:2] <- c("FID","IID") }
if(!(pheno_col %in% names(pheno))) quit(status=1)

data <- merge(geno[, .(FID, IID, high_risk, total_risk)], pheno[, c("FID","IID",pheno_col), with=FALSE], by=c("FID","IID"))
names(data)[names(data)==pheno_col] <- "PHENOTYPE"
data <- merge(data, covar, by=c("FID","IID"))

covar_list <- unlist(strsplit(covar_names, ","))
available_covars <- intersect(covar_list, names(data))
if(length(available_covars)>0) { formula_str <- paste("PHENOTYPE ~ high_risk +", paste(available_covars, collapse=" + ")) } else { formula_str <- "PHENOTYPE ~ high_risk" }

model <- tryCatch(lm(as.formula(formula_str), data=data), error=function(e) NULL)
if(!is.null(model) && "high_risk" %in% rownames(summary(model)$coefficients)) {
    coef <- summary(model)$coefficients["high_risk", ]
    results <- data.table(Test="APOL1_combined_recessive", Beta=coef["Estimate"], SE=coef["Std. Error"],
                          t_value=coef["t value"], P_value=coef["Pr(>|t|)"], N=nrow(model$model),
                          N_high_risk=sum(data$$high_risk==1, na.rm=TRUE), N_low_risk=sum(data$$high_risk==0, na.rm=TRUE))
    fwrite(results, file.path(results_dir, paste0(output_prefix, "_APOL1_combined_results.txt")), sep="\t")
}
EOF
            TOTAL_TESTS=$((TOTAL_TESTS + 1))
            echo "  Combined model: done"
        fi
    fi
else
    echo "  Combined model: skipped (missing G1 haplotype and G2)"
fi

# MYH9 tests
echo ""
echo "Running MYH9 tests..."

if [ "$MYH9_COUNT" -ge 1 ]; then
    plink2 --bfile $${GENO_PREFIX} --extract $${SNPLIST_DIR}/MYH9_snps.txt \
        --pheno $${PHENO_FILE} --pheno-name $${PHENO_COL} \
        --covar $${COVAR_FILE_WITH_PCS} --covar-name $${COVAR_NAMES} \
        --glm hide-covar --freq \
        --out $${RESULTS_DIR}/$${OUTPUT_PREFIX}_MYH9_additive
    TOTAL_TESTS=$((TOTAL_TESTS + MYH9_COUNT))
    echo "  Additive (${MYH9_COUNT} SNPs): done"
else
    echo "  Skipped (no SNPs available)"
fi

# Summary
if [ "$TOTAL_TESTS" -gt 0 ]; then
    BONF=$$(echo "scale=6; 0.05 / $${TOTAL_TESTS}" | bc)
else
    BONF="NA"
fi

cat >> ${SUMMARY_FILE} << EOF

Tests run: ${TOTAL_TESTS}
Bonferroni threshold: ${BONF}
EOF

rm -rf ${TEMP_DIR}

echo ""
echo "Complete. Results in ${RESULTS_DIR}"
echo "Tests: $${TOTAL_TESTS}, Bonferroni: $${BONF}"
