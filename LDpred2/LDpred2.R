# ==============================================================================
# LDpred2-auto workflow
# ==============================================================================
#
# Purpose
# -------
# This script constructs a polygenic score using LDpred2-auto, calculates scores
# in a target genotype dataset, and evaluates the incremental variance explained
# by the score.
#
# Method and software
# -------------------
# LDpred2 is implemented in the R package `bigsnpr`.
#
# Official repository:
# https://github.com/privefl/bigsnpr
#
# Official LDpred2 tutorial:
# https://privefl.github.io/bigsnpr/articles/LDpred2.html
#
# Package documentation:
# https://privefl.github.io/bigsnpr/
#
# Required local files
# --------------------
# 1. GWAS summary statistics:
#      path/to/project/sumstats/<study>.tsv
#
#    Required columns:
#      chr       chromosome
#      pos       base-pair position
#      rsid      variant identifier
#      a0        non-effect allele
#      a1        effect allele
#      beta      GWAS effect estimate
#      beta_se   standard error of beta
#      N         GWAS sample size or effective sample size
#
# 2. Phenotype and covariate file:
#      path/to/project/Phenofiles/<trait>_covars.tsv
#
#    Required columns:
#      IID, Phenotype, Age, Sex, PC1, ..., PC10
#
# 3. Target/reference genotype data converted to a bigSNP object:
#      path/to/project/LDpred2/LIFE_bed.rds
#      path/to/project/LDpred2/LIFE_bed.bk
#
#    Create these files once from PLINK BED/BIM/FAM files with:
#
#      bigsnpr::snp_readBed("path/to/genotypes.bed",
#                           backingfile = "path/to/project/LDpred2/LIFE_bed")
#
#    The `.rds` and `.bk` files must remain together in the same directory.
#
# 4. Genetic positions in centimorgans:
#      path/to/project/LDpred2/POS2.rds
#
#    This vector must correspond exactly to the variants and ordering in the
#    bigSNP genotype object. It can be generated with:
#
#      POS2 <- bigsnpr::snp_asGeneticPos(
#        CHR,
#        POS,
#        dir = file.path(project_dir, "LDpred2", "genetic_maps"),
#        ncores = n_cores
#      )
#      saveRDS(POS2, file.path(project_dir, "LDpred2", "POS2.rds"))
#
# Important
# ---------
# - All genotype and GWAS files must use compatible genome builds.
# - Summary statistics must be harmonized before running this script.
# - The LD reference should be ancestry-matched to the intended application.
# - Individual-level LIFE-Adult data are controlled-access and are not included
#   in the GitHub repository.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

library(bigsnpr)
library(bigstatsr)
library(data.table)
library(dplyr)
library(ggplot2)


# ------------------------------------------------------------------------------
# 2. User settings
# ------------------------------------------------------------------------------

# Change only this path to the root directory containing all project subfolders.
project_dir <- "path/to/project"

trait <- "Diabetes"
study <- "Mahajan.NatGenet2018b.T2D.European"

# Number of CPU cores used by bigsnpr.
n_cores <- 8

# Expected phenotype type:
#   "continuous" for quantitative traits
#   "binary" for case-control traits
phenotype_type <- "continuous"

# File paths derived from the single project prefix.
ldpred_dir <- file.path(project_dir, "LDpred2")
sumstats_file <- file.path(project_dir, "sumstats", paste0(study, ".tsv"))
phenotype_file <- file.path(
  project_dir,
  "Phenofiles",
  paste0(trait, "_covars.tsv")
)
bigsnp_file <- file.path(ldpred_dir, "LIFE_bed.rds")
genetic_pos_file <- file.path(ldpred_dir, "POS2.rds")
out_dir <- file.path(ldpred_dir, study)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 3. Check input files
# ------------------------------------------------------------------------------

required_files <- c(
  sumstats_file,
  phenotype_file,
  bigsnp_file,
  genetic_pos_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "The following required files were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}


# ------------------------------------------------------------------------------
# 4. Read phenotype data and GWAS summary statistics
# ------------------------------------------------------------------------------

phenotype_data <- fread(phenotype_file)
sumstats <- fread(sumstats_file)

required_sumstats_columns <- c(
  "chr", "pos", "rsid", "a0", "a1", "beta", "beta_se", "N"
)

missing_sumstats_columns <- setdiff(
  required_sumstats_columns,
  names(sumstats)
)

if (length(missing_sumstats_columns) > 0) {
  stop(
    "The summary-statistics file is missing these columns: ",
    paste(missing_sumstats_columns, collapse = ", ")
  )
}

required_phenotype_columns <- c(
  "IID", "Phenotype", "Age", "Sex", paste0("PC", 1:10)
)

missing_phenotype_columns <- setdiff(
  required_phenotype_columns,
  names(phenotype_data)
)

if (length(missing_phenotype_columns) > 0) {
  stop(
    "The phenotype file is missing these columns: ",
    paste(missing_phenotype_columns, collapse = ", ")
  )
}

n_sumstats_raw <- nrow(sumstats)


# ------------------------------------------------------------------------------
# 5. Attach the bigSNP genotype object
# ------------------------------------------------------------------------------

obj_bigSNP <- snp_attach(bigsnp_file)

G <- obj_bigSNP$genotypes
CHR <- obj_bigSNP$map$chromosome
POS <- obj_bigSNP$map$physical.pos

if (length(CHR) != ncol(G) || length(POS) != ncol(G)) {
  stop("The genotype map and genotype matrix have inconsistent dimensions.")
}


# ------------------------------------------------------------------------------
# 6. Match GWAS variants to the genotype data
# ------------------------------------------------------------------------------

# LDpred2 expects an effective sample-size column named `n_eff`.
sumstats[, n_eff := N]

# bigsnpr's map uses the allele order a1/a0 expected by snp_match().
genotype_map <- setNames(
  obj_bigSNP$map[, -3],
  c("chr", "rsid", "pos", "a1", "a0")
)

df_beta <- snp_match(sumstats, genotype_map)
n_df_beta_matched <- nrow(df_beta)

if (n_df_beta_matched == 0) {
  stop(
    "No variants matched between the summary statistics and genotype data. ",
    "Check genome builds, chromosome formatting, positions, and alleles."
  )
}

fwrite(
  df_beta,
  file = file.path(out_dir, "df_beta_matched_before_MAF.tsv.gz"),
  sep = "\t"
)


# ------------------------------------------------------------------------------
# 7. Load genetic positions
# ------------------------------------------------------------------------------

POS2 <- readRDS(genetic_pos_file)

if (length(POS2) != ncol(G)) {
  stop(
    "`POS2.rds` must contain one genetic position for every genotype variant ",
    "in the same order as the bigSNP object."
  )
}


# ------------------------------------------------------------------------------
# 8. Filter matched variants by minor allele frequency
# ------------------------------------------------------------------------------

ind_row <- rows_along(G)

maf <- snp_MAF(
  G,
  ind.row = ind_row,
  ind.col = df_beta$`_NUM_ID_`,
  ncores = n_cores
)

maf_threshold <- 1 / sqrt(length(ind_row))

df_beta[, maf := maf]
df_beta <- df_beta[maf > maf_threshold]

n_df_beta_after_maf <- nrow(df_beta)

if (n_df_beta_after_maf == 0) {
  stop("No variants remained after MAF filtering.")
}

fwrite(
  df_beta,
  file = file.path(out_dir, "df_beta_after_MAF.tsv.gz"),
  sep = "\t"
)

saveRDS(
  df_beta,
  file = file.path(out_dir, "df_beta_after_MAF.rds")
)


# ------------------------------------------------------------------------------
# 9. Create the sparse genome-wide LD correlation matrix
# ------------------------------------------------------------------------------

# LDpred2 should be run genome-wide. The correlation matrix is built chromosome
# by chromosome and combined into one on-disk sparse file-backed matrix.
sfbm_prefix <- tempfile(pattern = "ld_matrix_", tmpdir = out_dir)

corr <- NULL
ld <- numeric(0)

for (chr in 1:22) {
  
  message("Calculating LD for chromosome ", chr, "...")
  
  ind_chr <- which(df_beta$chr == chr)
  
  if (length(ind_chr) == 0) {
    message("No matched variants on chromosome ", chr, "; skipping.")
    next
  }
  
  genotype_indices <- df_beta$`_NUM_ID_`[ind_chr]
  
  corr_chr <- snp_cor(
    G,
    ind.col = genotype_indices,
    size = 3 / 1000,
    infos.pos = POS2[genotype_indices],
    ncores = n_cores
  )
  
  ld <- c(ld, Matrix::colSums(corr_chr^2))
  
  if (is.null(corr)) {
    corr <- as_SFBM(corr_chr, sfbm_prefix, compact = TRUE)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}

if (is.null(corr) || length(ld) == 0) {
  stop("The LD correlation matrix could not be created.")
}


# ------------------------------------------------------------------------------
# 10. Estimate SNP heritability using LD Score regression
# ------------------------------------------------------------------------------

ldsc <- with(
  df_beta,
  snp_ldsc(
    ld,
    length(ld),
    chi2 = (beta / beta_se)^2,
    sample_size = n_eff,
    blocks = NULL
  )
)

ldsc_h2_est <- ldsc[["h2"]]

if (!is.finite(ldsc_h2_est) || ldsc_h2_est <= 0) {
  stop(
    "LD Score regression returned an invalid heritability estimate: ",
    ldsc_h2_est
  )
}

saveRDS(
  ldsc,
  file = file.path(out_dir, "ldsc_results.rds")
)

write.table(
  as.data.frame(t(ldsc)),
  file = file.path(out_dir, "ldsc_results.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


# ------------------------------------------------------------------------------
# 11. Run LDpred2-auto
# ------------------------------------------------------------------------------

set.seed(1)

multi_auto <- snp_ldpred2_auto(
  corr = corr,
  df_beta = df_beta,
  h2_init = ldsc_h2_est,
  vec_p_init = seq_log(1e-4, 0.2, length.out = 30),
  allow_jump_sign = FALSE,
  shrink_corr = 0.4,
  use_MLE = FALSE,
  ncores = n_cores
)

saveRDS(
  multi_auto,
  file = file.path(out_dir, "multi_auto.rds")
)

chain_summary <- data.frame(
  chain = seq_along(multi_auto),
  h2_est = vapply(multi_auto, function(x) x$h2_est, numeric(1)),
  p_est = vapply(multi_auto, function(x) x$p_est, numeric(1)),
  corr_min = vapply(multi_auto, function(x) min(x$corr_est), numeric(1)),
  corr_max = vapply(multi_auto, function(x) max(x$corr_est), numeric(1))
)

fwrite(
  chain_summary,
  file = file.path(out_dir, "ldpred2_auto_chain_summary.tsv"),
  sep = "\t"
)


# ------------------------------------------------------------------------------
# 12. Save convergence plots for the first chain
# ------------------------------------------------------------------------------

auto1 <- multi_auto[[1]]

p_path_data <- data.frame(
  iteration = seq_along(auto1$path_p_est),
  p = auto1$path_p_est
)

h2_path_data <- data.frame(
  iteration = seq_along(auto1$path_h2_est),
  h2 = auto1$path_h2_est
)

p_path_plot <- ggplot(p_path_data, aes(x = iteration, y = p)) +
  geom_line() +
  theme_bigstatsr() +
  geom_hline(yintercept = auto1$p_est) +
  scale_y_log10() +
  labs(x = "Iteration", y = "Estimated proportion of causal variants")

h2_path_plot <- ggplot(h2_path_data, aes(x = iteration, y = h2)) +
  geom_line() +
  theme_bigstatsr() +
  geom_hline(yintercept = auto1$h2_est) +
  labs(x = "Iteration", y = "Estimated SNP heritability")

ggsave(
  file.path(out_dir, "chain1_p_path.png"),
  p_path_plot,
  width = 7,
  height = 4,
  dpi = 300
)

ggsave(
  file.path(out_dir, "chain1_h2_path.png"),
  h2_path_plot,
  width = 7,
  height = 4,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 13. Identify well-mixing chains
# ------------------------------------------------------------------------------

range_corr <- vapply(
  multi_auto,
  function(auto) diff(range(auto$corr_est)),
  numeric(1)
)

keep <- which(
  range_corr > 0.95 * quantile(range_corr, 0.95, na.rm = TRUE)
)

if (length(keep) == 0) {
  stop(
    "No LDpred2-auto chains passed the chain-filtering criterion. ",
    "Inspect the chain diagnostics before changing the threshold."
  )
}

fwrite(
  data.frame(
    chain = seq_along(range_corr),
    range_corr = range_corr,
    keep = seq_along(range_corr) %in% keep
  ),
  file = file.path(out_dir, "chain_ranges_and_keep.tsv"),
  sep = "\t"
)


# ------------------------------------------------------------------------------
# 14. Calculate final LDpred2 weights and individual scores
# ------------------------------------------------------------------------------

beta_matrix <- sapply(
  multi_auto[keep],
  function(auto) auto$beta_est
)

if (is.null(dim(beta_matrix))) {
  beta_auto <- as.numeric(beta_matrix)
} else {
  beta_auto <- rowMeans(beta_matrix)
}

pred_auto <- big_prodVec(
  G,
  beta_auto,
  ind.col = df_beta[["_NUM_ID_"]]
)

n_snps_used <- sum(beta_auto != 0)

message(
  "Number of variants with non-zero LDpred2 weights: ",
  n_snps_used
)

weights_table <- data.table(
  rsid = df_beta$rsid,
  chr = df_beta$chr,
  pos = df_beta$pos,
  a1 = df_beta$a1,
  a0 = df_beta$a0,
  maf = df_beta$maf,
  beta_auto = beta_auto
)

weights_prefix <- paste0("LDpred2_", trait, "_weights_beta_auto")

fwrite(
  weights_table,
  file = file.path(out_dir, paste0(weights_prefix, ".tsv.gz")),
  sep = "\t"
)

saveRDS(
  weights_table,
  file = file.path(out_dir, paste0(weights_prefix, ".rds"))
)


# ------------------------------------------------------------------------------
# 15. Merge the score with phenotype and covariate data
# ------------------------------------------------------------------------------

pred_auto_standardized <- as.numeric(scale(pred_auto))

prs_table <- data.frame(
  IID = obj_bigSNP$fam$sample.ID,
  PGS_LDpred2 = pred_auto_standardized
)

analysis_data <- phenotype_data |>
  left_join(prs_table, by = "IID") |>
  filter(
    complete.cases(
      Phenotype,
      Age,
      Sex,
      PC1,
      PC2,
      PC3,
      PC4,
      PC5,
      PC6,
      PC7,
      PC8,
      PC9,
      PC10,
      PGS_LDpred2
    )
  )

if (nrow(analysis_data) == 0) {
  stop("No complete samples remained after merging scores and covariates.")
}


# ------------------------------------------------------------------------------
# 16. Fit the prediction models
# ------------------------------------------------------------------------------

full_formula <- Phenotype ~ PGS_LDpred2 + Age + Sex +
  PC1 + PC2 + PC3 + PC4 + PC5 +
  PC6 + PC7 + PC8 + PC9 + PC10

reduced_formula <- Phenotype ~ Age + Sex +
  PC1 + PC2 + PC3 + PC4 + PC5 +
  PC6 + PC7 + PC8 + PC9 + PC10

if (phenotype_type == "continuous") {
  
  full_model <- lm(full_formula, data = analysis_data)
  reduced_model <- lm(reduced_formula, data = analysis_data)
  
  r2_full <- summary(full_model)$r.squared
  r2_covariates <- summary(reduced_model)$r.squared
  r2_incremental <- r2_full - r2_covariates
  
  performance_text <- paste0(
    "Covariate-only model R²: ", round(r2_covariates, 4), "\n",
    "Full model R²:           ", round(r2_full, 4), "\n",
    "Incremental PRS R²:      ", round(r2_incremental, 4), "\n"
  )
  
} else if (phenotype_type == "binary") {
  
  full_model <- glm(
    full_formula,
    data = analysis_data,
    family = binomial()
  )
  
  reduced_model <- glm(
    reduced_formula,
    data = analysis_data,
    family = binomial()
  )
  
  likelihood_ratio_test <- anova(
    reduced_model,
    full_model,
    test = "LRT"
  )
  
  performance_text <- paste0(
    "Binary phenotype: logistic regression was used.\n",
    "See likelihood_ratio_test.tsv for the nested-model comparison.\n"
  )
  
  fwrite(
    as.data.table(likelihood_ratio_test, keep.rownames = "model"),
    file = file.path(out_dir, "likelihood_ratio_test.tsv"),
    sep = "\t"
  )
  
} else {
  stop("`phenotype_type` must be either 'continuous' or 'binary'.")
}

fwrite(
  as.data.table(analysis_data),
  file = file.path(out_dir, paste0("LDpred2_", trait, "_scores.tsv")),
  sep = "\t"
)

capture.output(
  summary(full_model),
  file = file.path(out_dir, "full_model_summary.txt")
)

capture.output(
  summary(reduced_model),
  file = file.path(out_dir, "reduced_model_summary.txt")
)

cat(performance_text)
writeLines(
  performance_text,
  con = file.path(out_dir, "model_performance.txt")
)


# ------------------------------------------------------------------------------
# 17. Save run summary and session information
# ------------------------------------------------------------------------------

run_summary <- c(
  paste("Trait:", trait),
  paste("Study:", study),
  paste("Phenotype type:", phenotype_type),
  paste("Number of CPU cores:", n_cores),
  paste("N summary-statistic variants before matching:", n_sumstats_raw),
  paste("N matched variants:", n_df_beta_matched),
  paste("MAF threshold:", maf_threshold),
  paste("N variants after MAF filtering:", n_df_beta_after_maf),
  paste("N retained LDpred2-auto chains:", length(keep)),
  paste("N variants with non-zero final weights:", n_snps_used),
  paste("N samples in the final phenotype analysis:", nrow(analysis_data))
)

writeLines(
  run_summary,
  con = file.path(out_dir, "run_summary.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)