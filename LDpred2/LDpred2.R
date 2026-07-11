## ------------------------------------------------------------------
## Libraries
## ------------------------------------------------------------------
.libPaths("/net/ifs2/san_projekte/projekte/Doktoranden/MarieFreudenberg/forostar_Rlibs/")

# for forostarn number of BLAS-Threads needs to be set to 1 
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
library(RhpcBLASctl)
RhpcBLASctl::blas_set_num_threads(1)

library(bigsnpr)
library(bigstatsr)
library(data.table)
library(magrittr)
library(ggplot2)
library(readxl)
library(readr)
library(dplyr)
library(glue)

## ------------------------------------------------------------------
## Paths / trait
## ------------------------------------------------------------------
trait <- "Diabetes"
study <- "Mahajan.NatGenet2018b.T2D.European"
forostar <- TRUE

base_path <- if (forostar) {
  "/net/ifs2/san_projekte/projekte/Doktoranden/MarieFreudenberg/analysen/"
} else {
  "R:/Doktoranden/MarieFreudenberg/analysen/"
}

ldpred_dir <- file.path(base_path, "LDpred2")
out_dir    <- file.path(glue(base_path, "LDpred2/{study}"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------------
## Read phenotype + sumstats (sumstats need chr, pos, a0, a1, beta, beta_se, N)
## ------------------------------------------------------------------
data     <- fread(file.path(glue(base_path, "Phenofiles/{trait}_covars.tsv")))
sumstats <- fread(file.path(glue(base_path, "sumstats/{study}.tsv")))

n_sumstats_raw <- nrow(sumstats)

## ------------------------------------------------------------------
## Attach genotype (bigSNP)
## ------------------------------------------------------------------
obj.bigSNP <- snp_attach(file.path(ldpred_dir, "LIFE_bed.rds"))

G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection

NCORES <- 32  

## ------------------------------------------------------------------
## Match variants between sumstats and genotypes
## ------------------------------------------------------------------
sumstats$n_eff <- sumstats$N
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))

df_beta <- snp_match(sumstats, map)
df_beta_matched <- df_beta

n_df_beta_matched <- nrow(df_beta_matched)

# Save matched sumstats before MAF filter (optional but helpful)
fwrite(
  df_beta_matched,
  file = file.path(out_dir, "df_beta_matched_before_MAF.tsv.gz"),
  sep  = "\t"
)

## ------------------------------------------------------------------
## Genetic positions (cM) – reused from RDS
## ------------------------------------------------------------------
POS2 <- readRDS(file.path(ldpred_dir, "POS2.rds"))

## ------------------------------------------------------------------
## Filter on MAF in reference panel
## ------------------------------------------------------------------
ind.row <- rows_along(G)
maf     <- snp_MAF(G, ind.row = ind.row, ind.col = df_beta$`_NUM_ID_`, ncores = NCORES)
maf_thr <- 1 / sqrt(length(ind.row))

df_beta$maf <- maf
df_beta     <- df_beta[maf > maf_thr, ]

n_df_beta_maf <- nrow(df_beta)

# Save df_beta after MAF filtering
fwrite(
  df_beta,
  file = file.path(out_dir, "df_beta_after_MAF.tsv.gz"),
  sep  = "\t"
)
saveRDS(df_beta, file.path(out_dir, "df_beta_after_MAF.rds"))

## ------------------------------------------------------------------
## Create on-disk sparse genome-wide correlation matrix 
## ------------------------------------------------------------------

tmp <- tempfile(tmpdir = out_dir)

for (chr in 1:22) {
  
  print(chr)
  
  ind.chr <- which(df_beta$chr == chr)
  
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  
  corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000,
                   infos.pos = POS2[ind.chr2], ncores = NCORES)
  
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp, compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}


## ------------------------------------------------------------------
## LD Score regression
## ------------------------------------------------------------------
ldsc <- with(
  df_beta,
  snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
           sample_size = n_eff, blocks = NULL)
)

ldsc_h2_est <- ldsc[["h2"]]

saveRDS(ldsc, file.path(out_dir, "ldsc_results.rds"))
write.table(as.data.frame(t(ldsc)),
            file = file.path(out_dir, "ldsc_results.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

## ------------------------------------------------------------------
## LDpred2-auto
## ------------------------------------------------------------------
coef_shrink <- 0.95

set.seed(1)

multi_auto <- snp_ldpred2_auto(
  corr, df_beta,
  h2_init = ldsc_h2_est,
  vec_p_init = seq_log(1e-4, 0.2, length.out = 30),
  allow_jump_sign = FALSE,
  shrink_corr = 0.4,
  use_MLE = FALSE,
  ncores = NCORES
)

saveRDS(multi_auto, file.path(out_dir, "multi_auto.rds"))

chain_summary <- data.frame(
  chain    = seq_along(multi_auto),
  h2_est   = sapply(multi_auto, function(a) a$h2_est),
  p_est    = sapply(multi_auto, function(a) a$p_est),
  corr_min = sapply(multi_auto, function(a) min(a$corr_est)),
  corr_max = sapply(multi_auto, function(a) max(a$corr_est))
)
fwrite(chain_summary,
       file = file.path(out_dir, "ldpred2_auto_chain_summary.tsv"),
       sep  = "\t")

## ------------------------------------------------------------------
## Convergence plots for chain 1 (save PNGs)
## ------------------------------------------------------------------
auto1 <- multi_auto[[1]]

df_p  <- data.frame(iter = seq_along(auto1$path_p_est),
                    p    = auto1$path_p_est)
df_h2 <- data.frame(iter = seq_along(auto1$path_h2_est),
                    h2   = auto1$path_h2_est)

p_path <- ggplot(df_p, aes(x = iter, y = p)) +
  geom_line() +
  theme_bigstatsr() +
  geom_hline(yintercept = auto1$p_est, colour = "blue") +
  scale_y_log10() +
  labs(x = "Iteration", y = "p")

h2_path <- ggplot(df_h2, aes(x = iter, y = h2)) +
  geom_line() +
  theme_bigstatsr() +
  geom_hline(yintercept = auto1$h2_est, colour = "blue") +
  labs(x = "Iteration", y = "h2")

ggsave(file.path(out_dir, "chain1_p_path.png"),  p_path,  width = 7, height = 4, dpi = 300)
ggsave(file.path(out_dir, "chain1_h2_path.png"), h2_path, width = 7, height = 4, dpi = 300)

## ------------------------------------------------------------------
## Filter bad chains
## ------------------------------------------------------------------
range_corr <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
keep       <- which(range_corr > (0.95 * quantile(range_corr, 0.95, na.rm = TRUE)))

fwrite(
  data.frame(chain = seq_along(range_corr),
             range_corr = range_corr,
             keep = seq_along(range_corr) %in% keep),
  file = file.path(out_dir, "chain_ranges_and_keep.tsv"),
  sep  = "\t"
)

## ------------------------------------------------------------------
## PRS weights and scores
## ------------------------------------------------------------------
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
pred_auto <- big_prodVec(G, beta_auto, ind.col = df_beta[["_NUM_ID_"]])

num_snps_used <- sum(beta_auto != 0)
cat("Number of SNPs used in PRS prediction (non-zero weights):",
    num_snps_used, "\n")

weights_tbl <- data.table(
  rsid      = df_beta$rsid,
  chr       = df_beta$chr,
  pos       = df_beta$pos,
  a1        = df_beta$a1,
  a0        = df_beta$a0,
  maf       = df_beta$maf,
  beta_auto = beta_auto
)

fwrite(
  weights_tbl,
  file = file.path(out_dir, "LDpred2_BMI_weights_beta_auto.tsv.gz"),
  sep  = "\t"
)
saveRDS(
  weights_tbl,
  file.path(out_dir, "LDpred2_BMI_weights_beta_auto.rds")
)

## ------------------------------------------------------------------
## Attach scores to phenotype data, fit model
## ------------------------------------------------------------------
pred_auto_std <- as.numeric(scale(pred_auto))

prs <- tibble(
  IID = obj.bigSNP$fam$sample.ID,
  PGS_LDpred2 = pred_auto_std
)

data <- data %>%
  left_join(prs, by = "IID") %>%
  filter(complete.cases(Phenotype, Age, Sex,
                        PC1, PC2, PC3, PC4, PC5,
                        PC6, PC7, PC8, PC9, PC10,
                        PGS_LDpred2))

model <- lm(
  Phenotype ~ PGS_LDpred2 + Age + Sex +
    PC1 + PC2 + PC3 + PC4 + PC5 +
    PC6 + PC7 + PC8 + PC9 + PC10,
  data = data
)

write_tsv(
  data,
  file.path(out_dir, "LDPred2_BMI_scores_0.2h2.tsv")
)

capture.output(
  summary(model),
  file = file.path(out_dir, "model_summary.txt")
)


## ------------------------------------------------------------------
## Variance explained (R² decomposition)
## ------------------------------------------------------------------
model_reduced <- lm(Phenotype ~ Age + Sex +
                      PC1 + PC2 + PC3 + PC4 + PC5 +
                      PC6 + PC7 + PC8 + PC9 + PC10, data = data)

r2_full    <- summary(model)$r.squared
r2_covars <- summary(model_reduced)$r.squared
r2_prs <- r2_full - r2_covars
  
  var_explained_txt <- glue("
Base model R²:                {round(r2_covars, 4)}
Full model R²:                {round(r2_full, 4)}
Incremental R²:               {round(r2_prs, 4)}
")

cat(var_explained_txt, "\n")

writeLines(
  var_explained_txt,
  con = file.path(out_dir, "variance_explained_R2.txt")
)

## ------------------------------------------------------------------
## Small run summary + session info
## ------------------------------------------------------------------
run_summary <- c(
  paste("N_sumstats_raw:",        n_sumstats_raw),
  paste("N_df_beta_matched:",     n_df_beta_matched),
  paste("MAF_threshold:",         maf_thr),
  paste("N_df_beta_after_MAF:",   n_df_beta_maf),
  paste("N_snps_used_beta_auto:", num_snps_used)
)
writeLines(run_summary, file.path(out_dir, "run_summary.txt"))

capture.output(sessionInfo(),
               file = file.path(out_dir, "sessionInfo.txt"))
