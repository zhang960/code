MR：
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(TwoSampleMR)
})

exposureFile <- "D:/mendelian/mendelian/exposure.F.csv"
outcomeFile  <- "D:/mendelian/mendelian/mendelian/MetaSepsis_TwosampleMR.txt"
outcomeName  <- "A-U-F-MetaSepsis"
setwd("D:/mendelian/mendelian")

exp_all <- readr::read_csv(exposureFile, show_col_types = FALSE)

need_cols <- c("SNP","beta.exposure","se.exposure","pval.exposure",
               "effect_allele.exposure","other_allele.exposure",
               "eaf.exposure","exposure","id.exposure")
miss <- setdiff(need_cols, names(exp_all))
if (length(miss)) stop("Missing columns in exposure file: ", paste(miss, collapse=", "))

exp_all <- exp_all %>%
  filter(!is.na(SNP), !is.na(beta.exposure), !is.na(se.exposure), !is.na(pval.exposure))

if (nrow(exp_all) == 0) stop("Exposure data is empty or all values are NA.")

outcome_dat <- tryCatch({
  TwoSampleMR::read_outcome_data(
    snps     = unique(exp_all$SNP),
    filename = outcomeFile, sep = "\t",
    snp_col  = "SNP",
    beta_col = "beta",
    se_col   = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    pval_col = "pval",
    eaf_col  = "eaf"
  )
}, error = function(e) {
  message("Outcome file has no eaf or column names do not match: ", e$message)
  TwoSampleMR::read_outcome_data(
    snps     = unique(exp_all$SNP),
    filename = outcomeFile, sep = "\t",
    snp_col  = "SNP",
    beta_col = "beta",
    se_col   = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    pval_col = "pval"
  )
})
outcome_dat$outcome <- outcomeName

groups <- split(exp_all, exp_all$exposure)

mr_rows  <- list()
snp_rows <- list()
k <- 0; m <- 0

for (ex_name in names(groups)) {
  ex_df <- groups[[ex_name]]
  message("\n==============================")
  message("Exposure: ", ex_name, "  raw SNPs=", nrow(ex_df))
  
  ex_df <- ex_df %>% arrange(pval.exposure) %>% distinct(SNP, .keep_all = TRUE)
  
  out_sub <- outcome_dat %>% filter(SNP %in% ex_df$SNP)
  if (nrow(out_sub) == 0) { message("No matched SNPs in outcome data, skipped"); next }
  
  dat <- tryCatch(TwoSampleMR::harmonise_data(exposure_dat = ex_df, outcome_dat = out_sub, action = 2),
                  error=function(e){ message("Harmonisation failed: ", e$message); NULL })
  if (is.null(dat) || nrow(dat)==0) { message("No data after harmonisation, skipped"); next }
  
  keep <- dat %>% filter(mr_keep == TRUE)
  if (nrow(keep) == 0) { message("No valid IVs, skipped"); next }
  
  m <- m + 1
  snp_rows[[m]] <- keep %>%
    mutate(exposure_name = ex_name) %>%
    select(exposure_name, SNP, effect_allele.exposure, other_allele.exposure,
           beta.exposure, se.exposure, eaf.exposure, pval.exposure)
  
  nsnp <- dplyr::n_distinct(keep$SNP)
  use_method <- "mr_ivw_fe"
  Q <- Q_pval <- NA_real_
  
  if (nsnp == 1) {
    message("Single SNP: Wald ratio")
    use_method <- "mr_wald_ratio"
  } else {
    het <- tryCatch(TwoSampleMR::mr_heterogeneity(keep, method_list = c("mr_ivw")),
                    error=function(e) NULL)
    if (!is.null(het) && nrow(het)>0) {
      Q <- het$Q[1]; Q_pval <- het$Q_pval[1]
      if (!is.na(Q_pval) && Q_pval < 0.05) {
        message(sprintf("Significant heterogeneity (Q=%.3f, p=%.3g): using IVW-MRE", Q, Q_pval))
        use_method <- "mr_ivw_mre"
      } else {
        message(sprintf("No significant heterogeneity (Q=%.3f, p=%.3g): using IVW-FE", Q, Q_pval))
        use_method <- "mr_ivw_fe"
      }
    } else {
      message("No heterogeneity result: using IVW-FE by default")
    }
  }
  
  mr_res <- tryCatch(TwoSampleMR::mr(keep, method_list = use_method),
                     error=function(e){ message("MR failed: ", e$message); NULL })
  if (is.null(mr_res) || nrow(mr_res)==0) { message("MR result is empty, skipped"); next }
  
  or_res <- tryCatch(TwoSampleMR::generate_odds_ratios(mr_res), error=function(e) NULL)
  if (!is.null(or_res)) {
    mr_res$OR    <- or_res$or
    mr_res$OR_lo <- or_res$or_lci95
    mr_res$OR_hi <- or_res$or_uci95
  } else {
    mr_res$OR <- mr_res$OR_lo <- mr_res$OR_hi <- NA_real_
  }
  
  k <- k + 1
  mr_rows[[k]] <- mr_res %>%
    mutate(
      exposure_name = ex_name,
      nsnp = nsnp,
      method_used = use_method,
      het_Q = Q, het_Q_pval = Q_pval
    ) %>%
    select(exposure_name, outcome, id.exposure, method, method_used,
           nsnp, b, se, pval, OR, OR_lo, OR_hi)
}

if (length(mr_rows) > 0) {
  MR_sum <- bind_rows(mr_rows)
  write_csv(MR_sum, "MR_by_exposure.csv")
  message("Saved: MR_by_exposure.csv")
}

if (length(snp_rows) > 0) {
  SNP_sum <- bind_rows(snp_rows)
  write_csv(SNP_sum, "SNPs_by_exposure.csv")
  message("Saved: SNPs_by_exposure.csv")
}
Reverse MR:
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(stringr)

setwd("D:\\data\\metabolites")   
exposure_file <- "exposure.F.csv"            
outcome_file  <- "outcomeID.txt"          

log_msg   <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse=" ")))
safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

map_to_rsids <- function(snp_vec) {
  snp_vec <- unique(snp_vec)
  is_rs   <- grepl("^rs\\d+$", snp_vec)
  keep_rs <- snp_vec[is_rs]
  to_map  <- snp_vec[!is_rs]
  if (!length(to_map)) return(unique(keep_rs))
  log_msg("Number of non-rsID variants to map: ", length(to_map))
  mv <- tryCatch(ieugwasr::map_variants(to_map), error=function(e) NULL)
  if (is.null(mv) || !"rsid" %in% names(mv)) {
    log_msg("Mapping failed or returned empty; only existing rsIDs are retained")
    return(unique(keep_rs))
  }
  mapped_rs <- unique(na.omit(mv$rsid))
  log_msg("Number of variants successfully mapped to rsID: ", length(mapped_rs))
  unique(c(keep_rs, mapped_rs))
}

extract_outcome_with_retry <- function(snps, outcome_id, proxies=TRUE, rsq=0.8,
                                       max_retries_other=3) {
  attempt <- 0
  repeat {
    attempt <- attempt + 1
    res <- tryCatch(
      TwoSampleMR::extract_outcome_data(
        snps          = snps,
        outcomes      = outcome_id,
        proxies       = proxies,
        rsq           = rsq,
        align_alleles = 0,
        palindromes   = 1,
        maf_threshold = 0.3
      ),
      error = function(e) e
    )
    if (!inherits(res, "error")) return(res)
    
    msg <- conditionMessage(res)
    
    if (grepl("used up your OpenGWAS allowance", msg, ignore.case = TRUE)) {
      ts_str <- tryCatch({
        m <- regexpr("\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}(?:\\.\\d+)?", msg, perl = TRUE)
        if (m[1] > 0) regmatches(msg, m)[[1]] else NA_character_
      }, error = function(e) NA_character_)
      t_target <- tryCatch({
        if (!is.na(ts_str)) as.POSIXct(ts_str, tz = Sys.timezone()) else NA
      }, warning = function(w) NA, error = function(e) NA)
      wait_sec <- if (!is.na(t_target)) {
        max(2, as.numeric(difftime(t_target, Sys.time(), units = "secs")))
      } else 90
      cat(sprintf("[429] OpenGWAS rate limit: wait ~%.0f seconds and retry (%s)\n",
                  wait_sec, if (!is.na(t_target)) format(t_target, "%H:%M:%S") else "no reset time parsed"))
      Sys.sleep(wait_sec + 2)
      next
    }
    
    if (attempt <= max_retries_other) {
      cat(sprintf("[Warning] Extraction failed (attempt %d/%d): %s; retrying after 10 seconds.\n",
                  attempt, max_retries_other, msg))
      Sys.sleep(10)
      next
    }
    
    cat(sprintf("[Error] Failed after multiple retries: %s; extraction skipped.\n", msg))
    return(NULL)
  }
}

run_mr_for_one_outcome_by_id <- function(exp_batch, out_dat_all, outcome_id) {
  group_list <- split(exp_batch, exp_batch$id.exposure)
  res_keep <- data.frame()
  
  gi <- 0
  for (gid in names(group_list)) {
    gi <- gi + 1
    exp_dat <- group_list[[gid]]
    exposure_id   <- unique(exp_dat$id.exposure)
    exposure_name <- unique(exp_dat$exposure)
    prefix <- paste0(safe_name(exposure_name), "__", safe_name(exposure_id), "_", safe_name(outcome_id))
    
    log_msg(sprintf("Group %d/%d: id.exposure=%s | exposure=%s | SNP=%d",
                    gi, length(group_list), exposure_id, exposure_name, length(unique(exp_dat$SNP))))
    
    out_dat <- out_dat_all[out_dat_all$SNP %in% exp_dat$SNP, ]
    if (!nrow(out_dat)) { 
      log_msg("No overlapping outcome data, skipped") 
      next 
    }
    
    dat <- TwoSampleMR::harmonise_data(exp_dat, out_dat, action = harmonise_action)
    
    if (!nrow(dat) || !any(dat$mr_keep)) {
      log_msg("No valid IVs after harmonisation, skipped")
      next
    }
    
    dat <- dat[dat$mr_keep, , drop = FALSE]
    n_iv <- nrow(dat)
    if (n_iv == 0) { 
      log_msg("IV count is 0 after filtering, skipped") 
      next 
    }
    
    mr_results <- tryCatch(TwoSampleMR::mr(dat), error = function(e) {
      log_msg("MR error: ", conditionMessage(e), " skipped")
      return(NULL)
    })
    if (is.null(mr_results) || !nrow(mr_results)) { 
      log_msg("MR returned empty result, skipped") 
      next 
    }
    
    mr_results <- mr_results[is.finite(mr_results$b) & is.finite(mr_results$se), , drop = FALSE]
    if (!nrow(mr_results)) { 
      log_msg("MR results are all NA or non-finite, skipped") 
      next 
    }
    
    mr_or <- tryCatch(TwoSampleMR::generate_odds_ratios(mr_results),
                      error = function(e) { 
                        log_msg("OR generation error: ", conditionMessage(e)); 
                        NULL 
                      })
    if (is.null(mr_or) || !nrow(mr_or)) { 
      log_msg("OR result is empty, skipped") 
      next 
    }
    
    mr_or$n_iv <- n_iv
    mr_or$outcome_id <- outcome_id
    
    res_keep <- dplyr::bind_rows(res_keep, mr_or)
    write.csv(mr_or, paste0(prefix, "_MR_results.csv"), row.names = FALSE)
    log_msg("Saved: ", prefix, "_MR_results.csv")
    
    outTab <- dat[dat$mr_keep == TRUE, ]
    write.csv(outTab, paste0(prefix, ".table.SNP.csv"), row.names = FALSE)
    log_msg("Saved: ", prefix, ".table.SNP.csv")
    
    if (n_iv > 1 && !"Wald ratio" %in% mr_results$method) {
      try(write.csv(TwoSampleMR::mr_heterogeneity(dat), paste0(prefix, "_heterogeneity.csv"), row.names = FALSE), silent = TRUE)
      try(write.csv(TwoSampleMR::mr_pleiotropy_test(dat), paste0(prefix, "_pleiotropy.csv"), row.names = FALSE), silent = TRUE)
      try({
        p1 <- TwoSampleMR::mr_scatter_plot(mr_results, dat)
        pdf(paste0(prefix, "_scatter_plot.pdf"), 7, 6.5); print(p1[[1]]); dev.off()
      }, silent = TRUE)
      try({
        res_single <- TwoSampleMR::mr_singlesnp(dat)
        p2 <- TwoSampleMR::mr_forest_plot(res_single)
        pdf(paste0(prefix, "_forest_plot.pdf"), 6.5, 5); print(p2[[1]]); dev.off()
      }, silent = TRUE)
      try({
        res_single <- TwoSampleMR::mr_singlesnp(dat)
        p3 <- TwoSampleMR::mr_funnel_plot(res_single)
        pdf(paste0(prefix, "_funnel_plot.pdf"), 6.5, 6); print(p3[[1]]); dev.off()
      }, silent = TRUE)
      try({
        p4 <- TwoSampleMR::mr_leaveoneout_plot(TwoSampleMR::mr_leaveoneout(dat))
        pdf(paste0(prefix, "_leaveoneout_plot.pdf"), 6.5, 5); print(p4[[1]]); dev.off()
      }, silent = TRUE)
      log_msg("Sensitivity plots exported")
    } else {
      log_msg("Single IV or insufficient IVs, sensitivity plots skipped")
    }
  }
  res_keep
}

log_msg("Reading exposure file: ", exposure_file)
exposure_data <- read.csv(exposure_file, header = TRUE)

log_msg("Reading outcome IDs: ", outcome_file)
outcome_ids <- readLines(outcome_file)
outcome_ids <- outcome_ids[nchar(trimws(outcome_ids)) > 0]
if (length(outcome_ids) == 0) stop("outcomeID.txt is empty")
log_msg("Number of outcomes: ", length(outcome_ids))

log_msg("Number of raw exposure SNPs: ", length(unique(exposure_data$SNP)))
mapped_rs <- map_to_rsids(exposure_data$SNP)

exposure_data$SNP <- ifelse(grepl("^rs\\d+$", exposure_data$SNP), exposure_data$SNP, NA_character_)
exposure_data <- exposure_data[exposure_data$SNP %in% mapped_rs, ]
log_msg("Number of exposure rows after mapping: ", nrow(exposure_data))

exp_all <- TwoSampleMR::format_data(
  exposure_data, type="exposure", snps=NULL, header=TRUE,
  phenotype_col="exposure", snp_col="SNP",
  beta_col="beta.exposure", se_col="se.exposure",
  eaf_col="eaf.exposure",
  effect_allele_col="effect_allele.exposure",
  other_allele_col="other_allele.exposure",
  pval_col="pval.exposure",
  samplesize_col="samplesize.exposure",
  id_col="id.exposure"
)

group_ids <- unique(exp_all$id.exposure)
log_msg("Total number of id.exposure: ", length(group_ids))
batch_idx <- ceiling(seq_along(group_ids) / groups_per_batch)
id_batches <- split(group_ids, batch_idx)
log_msg("Number of batches: ", length(id_batches), "; target id.exposure per batch: ", groups_per_batch)

all_combined <- data.frame()

for (b in seq_along(id_batches)) {
  ids_in_batch <- id_batches[[b]]
  exp_batch <- exp_all[exp_all$id.exposure %in% ids_in_batch, ]
  batch_snps <- unique(exp_batch$SNP)
  
  log_msg(sprintf("Start batch %d/%d: id.exposure=%d | SNP=%d",
                  b, length(id_batches), length(ids_in_batch), length(batch_snps)))
  
  for (oid in outcome_ids) {
    log_msg("Outcome: ", oid, " batch ", b)
    out_dat <- extract_outcome_with_retry(batch_snps, oid, proxies=proxies, rsq=proxy_r2)
    if (is.null(out_dat) || !nrow(out_dat)) {
      log_msg("No outcome data in this batch, skipped")
      next
    }
    res_one <- run_mr_for_one_outcome_by_id(exp_batch, out_dat, oid)
    if (!is.null(res_one) && nrow(res_one)) {
      res_one$batch <- b
      all_combined <- dplyr::bind_rows(all_combined, res_one)
      write.csv(all_combined, "combined_MR_results.csv", row.names = FALSE)
      log_msg("Saved cumulative summary: combined_MR_results.csv")
    }
    Sys.sleep(sleep_per_outcome)
  }
  
  log_msg(sprintf("End batch %d/%d", b, length(id_batches)))
  Sys.sleep(sleep_per_batch)
}

log_msg("All analyses completed. Summary file: combined_MR_results.csv")

GEO：
setwd("D:\\GEO\\")

library(limma)
library(impute)
library(openxlsx)

geo_raw <- read.table("geneMatrix.txt", sep = "\t", header = TRUE, check.names = FALSE)
geo_raw <- as.matrix(geo_raw)
rownames(geo_raw) <- geo_raw[, 1]
geo_exp <- geo_raw[, -1, drop = FALSE]

geo_exp <- apply(geo_exp, 2, as.numeric)
rownames(geo_exp) <- rownames(geo_raw)

geo_exp_imp <- impute.knn(geo_exp)$data
geo_data <- avereps(geo_exp_imp)

pdf("raw_box.pdf")
boxplot(geo_data, col = "green", xaxt = "n", outline = FALSE)
dev.off()

geo_data <- normalizeBetweenArrays(geo_data)

pdf("normal_box.pdf")
boxplot(geo_data, col = "red", xaxt = "n", outline = FALSE)
dev.off()

class <- c(rep("normal", 10), rep("treatment", 29))
stopifnot(length(class) == ncol(geo_data))

class_factor <- factor(class, levels = c("normal", "treatment"))
design <- model.matrix(~0 + class_factor)
colnames(design) <- levels(class_factor)

print(colSums(design))

fit <- lmFit(geo_data, design)
cont.matrix <- makeContrasts(treatment - normal, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allgene <- topTable(fit2, adjust.method = "fdr", number = Inf)

diffgene <- allgene[abs(allgene$logFC) >= 1 & allgene$adj.P.Val < 0.05, ]
Upgene   <- allgene[allgene$logFC >= 1 & allgene$adj.P.Val < 0.05, ]
Downgene <- allgene[allgene$logFC <= -1 & allgene$adj.P.Val < 0.05, ]

diffexp <- geo_data[rownames(diffgene), , drop = FALSE]

wb <- createWorkbook()

addWorksheet(wb, "allgene")
writeData(wb, "allgene", cbind(GeneID = rownames(allgene), allgene))

addWorksheet(wb, "diffgene")
writeData(wb, "diffgene", cbind(GeneID = rownames(diffgene), diffgene))

addWorksheet(wb, "upgene")
writeData(wb, "upgene", cbind(GeneID = rownames(Upgene), Upgene))

addWorksheet(wb, "down")
writeData(wb, "down", cbind(GeneID = rownames(Downgene), Downgene))

addWorksheet(wb, "diffexp")
writeData(wb, "diffexp", cbind(GeneID = rownames(diffexp), diffexp))

saveWorkbook(wb, "limma_results.xlsx", overwrite = TRUE)

cat("Done! Excel saved as limma_results.xlsx in: ", getwd(), "\n")


machine learning:
rm(list = ls())
set.seed(825)

library(readxl)
library(dplyr)
library(caret)
library(ranger)
library(xgboost)
library(SHAPforxgboost)
library(Boruta)
library(gbm)
library(kernlab)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(viridis)
library(tibble)

setwd("D:/ML")
data_raw <- read_excel("example_data(1).xlsx")

colnames(data_raw)[1:39] <- c(
  "id", "label", "Sex", "Age", "Weight", "BMI",
  "APACHE_II_24h", "SOFA", "GCS", "MAP",
  "Diabetic", "kidney_stone", "CHD", "Cerebrovascular",
  "ECMO", "Dementia", "Copd", "Paralysis", "CKD",
  "oXiris", "CRRT", "Mech_vent", "WBC", "NEU", "LYM",
  "ALB", "TBIL", "DBIL", "IBIL", "PLT", "PDW", "MPV",
  "VASO", "GLU", "HDL", "LDL", "TC", "TG", "TYG"
)

colnames(data_raw)[40] <- "Creat"

data <- data_raw

clean_numeric <- function(x){
  x <- as.character(x)
  x <- gsub(",", ".", x)
  x <- gsub("[^0-9\\.]", "", x)
  x[x == ""] <- NA
  as.numeric(x)
}

data$HDL <- clean_numeric(data$HDL)
data$LDL <- clean_numeric(data$LDL)
data$TC  <- clean_numeric(data$TC)
data$TG  <- clean_numeric(data$TG)
data$TYG <- clean_numeric(data$TYG)

data$label <- as.numeric(data$label)

data$label_factor <- factor(
  data$label,
  levels = c(0, 1),
  labels = c("control", "sepsis")
)

predictors <- setdiff(names(data), c("id", "label", "label_factor"))

nzv_info <- nearZeroVar(data[, predictors], saveMetrics = TRUE)
zero_var_cols <- rownames(nzv_info[nzv_info$zeroVar == TRUE, , drop = FALSE])
predictors <- setdiff(predictors, zero_var_cols)

preProc <- preProcess(
  data[, predictors],
  method = c("medianImpute", "center", "scale")
)

X_all <- predict(preProc, data[, predictors])
y_all <- data$label
y_fac <- data$label_factor

top_n <- 10

plot_importance_scaled <- function(df, var_col, imp_col, title, fill_color){
  df <- df %>%
    mutate(Scaled = .data[[imp_col]] / max(.data[[imp_col]], na.rm = TRUE))
  
  ggplot(df,
         aes(x = reorder(.data[[var_col]], Scaled),
             y = Scaled)) +
    geom_col(width = 0.7, fill = fill_color) +
    coord_flip() +
    theme_bw(base_size = 14) +
    theme(panel.grid.major.y = element_blank()) +
    labs(x = "Variable", y = "Relative importance (0-1)",
         title = title)
}

scale01 <- function(x) x / max(x, na.rm = TRUE)

set.seed(825)
bor <- Boruta(
  x = X_all,
  y = y_fac,
  doTrace = 0
)

bor_final <- TentativeRoughFix(bor)

bor_stats <- attStats(bor_final) %>%
  rownames_to_column("Variable") %>%
  dplyr::select(Variable, meanImp, medianImp, decision)

bor_imp <- bor_stats %>%
  filter(decision %in% c("Confirmed", "Tentative")) %>%
  arrange(desc(medianImp)) %>%
  mutate(Boruta = medianImp)

boruta_vars <- bor_imp$Variable

bor_top <- bor_imp[1:min(top_n, nrow(bor_imp)), c("Variable", "Boruta")]
colnames(bor_top)[1] <- "var"

imp_hist <- as.data.frame(bor_final$ImpHistory)
real_vars <- bor_stats$Variable
imp_real  <- imp_hist[, real_vars, drop = FALSE]

bor_long <- imp_real %>%
  pivot_longer(everything(),
               names_to = "Variable",
               values_to = "Importance") %>%
  left_join(bor_stats[, c("Variable", "decision")],
            by = "Variable") %>%
  na.omit()

bor_colors <- c(
  Confirmed = "#1b9e77",
  Tentative = "#7570b3",
  Rejected  = "#d95f02"
)

p_bor <- ggplot(bor_long,
                aes(x = reorder(Variable, Importance, FUN = median),
                    y = Importance,
                    fill = decision)) +
  geom_boxplot(
    color          = "black",
    alpha          = 0.95,
    outlier.size   = 0.8,
    outlier.alpha  = 0.55,
    outlier.stroke = 0.1
  ) +
  scale_fill_manual(values = bor_colors, na.translate = FALSE) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x        = element_text(angle = 60, hjust = 1, vjust = 1, size = 9),
    panel.grid.major.x = element_blank(),
    legend.position    = "right",
    plot.title         = element_text(size = 20, face = "bold")
  ) +
  labs(
    title = "Boruta Feature Selection",
    x     = "Clinical Variable",
    y     = "Importance across Boruta Iterations",
    fill  = "Decision"
  )

ggsave("Fig_Boruta_boxplot.pdf",
       p_bor, width = 30, height = 16, units = "cm", dpi = 600)

X_sel <- X_all[, boruta_vars, drop = FALSE]

set.seed(825)
rf_model <- ranger(
  y_fac ~ .,
  data       = data.frame(y_fac = y_fac, X_sel),
  num.trees  = 500,
  probability = FALSE,
  importance = "impurity"
)

rf_imp <- data.frame(
  var    = names(rf_model$variable.importance),
  RF_imp = as.numeric(rf_model$variable.importance)
) %>% arrange(desc(RF_imp))

rf_top <- rf_imp[1:min(top_n, nrow(rf_imp)), ]

p_rf <- plot_importance_scaled(
  rf_top, "var", "RF_imp",
  "Random Forest (Gini, Boruta-selected features)",
  "#4DBBD5FF"
)

ggsave("Fig_RF_importance.pdf",
       p_rf, width = 18, height = 12, units = "cm", dpi = 400)

set.seed(825)
rf_extratrees <- ranger(
  y_fac ~ .,
  data       = data.frame(y_fac = y_fac, X_sel),
  importance = "impurity",
  splitrule  = "extratrees",
  num.trees  = 500
)

et_imp <- data.frame(
  var    = names(rf_extratrees$variable.importance),
  ET_imp = as.numeric(rf_extratrees$variable.importance)
) %>% arrange(desc(ET_imp))

et_top <- et_imp[1:min(top_n, nrow(et_imp)), ]

p_et <- plot_importance_scaled(
  et_top, "var", "ET_imp",
  "Extremely Randomized Trees (splitrule = 'extratrees')",
  "#8491B4FF"
)

ggsave("Fig_ExtraTrees_importance.pdf",
       p_et, width = 18, height = 12, units = "cm", dpi = 400)

dtrain_xgb <- xgb.DMatrix(as.matrix(X_sel), label = y_all)

set.seed(825)
xgb_model <- xgb.train(
  params = list(
    objective        = "binary:logistic",
    eval_metric      = "auc",
    eta              = 0.05,
    max_depth        = 4,
    min_child_weight = 3,
    subsample        = 0.8,
    colsample_bytree = 0.8
  ),
  data    = dtrain_xgb,
  nrounds = 200,
  verbose = 0
)

xgb_shap_values <- shap.values(
  xgb_model = xgb_model,
  X_train   = as.matrix(X_sel)
)

xgb_shap_long <- shap.prep(
  shap_contrib = xgb_shap_values$shap_score,
  X_train      = as.matrix(X_sel)
)

xgb_imp <- data.frame(
  var      = names(xgb_shap_values$mean_shap_score),
  XGB_SHAP = as.numeric(xgb_shap_values$mean_shap_score)
) %>% arrange(desc(XGB_SHAP))

xgb_top <- xgb_imp[1:min(top_n, nrow(xgb_imp)), ]

p_xgb <- plot_importance_scaled(
  xgb_top, "var", "XGB_SHAP",
  "XGBoost-SHAP (mean |SHAP|, Boruta-selected features)",
  "#00A087FF"
)

ggsave("Fig_XGB_SHAP_bar.pdf",
       p_xgb, width = 18, height = 12, units = "cm", dpi = 400)

p_xgb_summary <- shap.plot.summary(xgb_shap_long)
ggsave("Fig_XGB_SHAP_summary.pdf",
       p_xgb_summary, width = 20, height = 16, units = "cm", dpi = 400)

gbm_data <- data.frame(y = y_all, X_sel)

set.seed(825)
gbm_model <- gbm(
  formula           = y ~ .,
  data              = gbm_data,
  distribution      = "bernoulli",
  n.trees           = 2000,
  interaction.depth = 3,
  shrinkage         = 0.01,
  bag.fraction      = 0.8,
  n.minobsinnode    = 5,
  train.fraction    = 1.0,
  verbose           = FALSE
)

gbm_imp_raw <- summary(gbm_model, plotit = FALSE)

gbm_imp <- gbm_imp_raw %>%
  as.data.frame() %>%
  rename(var = var, GBM_imp = rel.inf) %>%
  arrange(desc(GBM_imp))

gbm_top <- gbm_imp[1:min(top_n, nrow(gbm_imp)), ]

p_gbm <- plot_importance_scaled(
  gbm_top, "var", "GBM_imp",
  "GBM (Gradient Boosting) Importance",
  "#F39B7FFF"
)

ggsave("Fig_GBM_importance.pdf",
       p_gbm, width = 18, height = 12, units = "cm", dpi = 400)

set.seed(825)
ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

svm_fit <- train(
  y_fac ~ .,
  data = data.frame(y_fac = y_fac, X_sel),
  method = "svmLinear",
  trControl = ctrl,
  metric   = "ROC"
)

svm_imp_raw <- varImp(svm_fit, scale = FALSE)$importance
svm_col <- names(svm_imp_raw)[1]

svm_imp <- svm_imp_raw %>%
  rownames_to_column("var") %>%
  dplyr::rename(SVM_imp = !!svm_col) %>%
  arrange(desc(SVM_imp))

svm_top <- svm_imp[1:min(top_n, nrow(svm_imp)), ]

p_svm <- plot_importance_scaled(
  svm_top, "var", "SVM_imp",
  "Linear SVM (svmLinear) Importance",
  "#3C5488FF"
)

ggsave("Fig_SVM_importance.pdf",
       p_svm, width = 18, height = 12, units = "cm", dpi = 400)

vars_bor  <- bor_top$var
vars_rf   <- rf_top$var
vars_et   <- et_top$var
vars_xgb  <- xgb_top$var
vars_gbm  <- gbm_top$var
vars_svm  <- svm_top$var

all_list <- list(
  Boruta     = vars_bor,
  RF         = vars_rf,
  ExtraTrees = vars_et,
  XGB        = vars_xgb,
  GBM        = vars_gbm,
  SVM        = vars_svm
)

vars_intersect_6 <- Reduce(intersect, all_list)
vars_intersect_6

all_vars <- unique(unlist(all_list))
count_models <- sapply(all_vars, function(v){
  sum(sapply(all_list, function(vec) v %in% vec))
})

vars_in_3plus <- names(count_models[count_models >= 3])
vars_in_3plus_boruta <- intersect(vars_in_3plus, boruta_vars)
vars_in_3plus_boruta

imp_all <- bor_imp %>%
  dplyr::select(var = Variable, Boruta) %>%
  full_join(rf_imp,    by = "var") %>%
  full_join(et_imp,    by = "var") %>%
  full_join(xgb_imp,   by = "var") %>%
  full_join(gbm_imp,   by = "var") %>%
  full_join(svm_imp %>% dplyr::select(var, SVM_imp), by = "var")

imp_all_scaled <- imp_all %>%
  mutate(
    Boruta_s = scale01(Boruta),
    RF_s     = scale01(RF_imp),
    ET_s     = scale01(ET_imp),
    XGB_s    = scale01(XGB_SHAP),
    GBM_s    = scale01(GBM_imp),
    SVM_s    = scale01(SVM_imp)
  ) %>%
  dplyr::select(var, Boruta_s, RF_s, ET_s, XGB_s, GBM_s, SVM_s)

imp_long <- imp_all_scaled %>%
  pivot_longer(cols = c("Boruta_s", "RF_s", "ET_s", "XGB_s", "GBM_s", "SVM_s"),
               names_to = "Model",
               values_to = "Scaled") %>%
  filter(var %in% vars_in_3plus_boruta)

p_heat <- ggplot(imp_long,
                 aes(x = Model, y = var, fill = Scaled)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0.5
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Feature Importance Across 6 Methods",
    x = "",
    y = "",
    fill = "Scaled importance"
  )

ggsave("Fig_MultiModel_importance_heatmap_6methods_SVM.pdf",
       p_heat, width = 22, height = 18, units = "cm", dpi = 400)

p_panel_7 <- grid.arrange(
  p_bor,
  p_rf,
  p_et,
  p_xgb,
  p_gbm,
  p_svm,
  p_heat,
  ncol = 2
)

ggsave("Fig_Boruta_RF_ET_XGB_GBM_SVM_Heatmap_panel.pdf",
       p_panel_7, width = 32, height = 24, units = "cm", dpi = 400)

presence_df <- data.frame(
  var = rep(all_vars, each = length(all_list)),
  Model = rep(names(all_list), times = length(all_vars)),
  selected = 0
)

for(m in names(all_list)){
  presence_df$selected[
    presence_df$Model == m & presence_df$var %in% all_list[[m]]
  ] <- 1
}

core_vars <- vars_intersect_6

core_df <- presence_df %>%
  filter(var %in% core_vars)

core_df$Model <- factor(core_df$Model,
                        levels = c("Boruta", "RF", "ExtraTrees", "XGB", "GBM", "SVM"))

p_core_bar <- ggplot(core_df,
                     aes(x = Model, y = selected, fill = Model)) +
  geom_col(width = 0.7) +
  facet_wrap(~ var, ncol = 2) +
  scale_y_continuous(breaks = 0:1, limits = c(0, 1.05)) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    strip.background   = element_rect(fill = "grey95"),
    strip.text         = element_text(face = "bold")
  ) +
  labs(
    title = "Core variables selected by 6 ML models (0/1)",
    x     = "Model",
    y     = "Selected (0/1)"
  )

ggsave("Fig_Core4_modelcolor.pdf",
       p_core_bar, width = 20, height = 14, units = "cm", dpi = 400)

layout_mat <- matrix(
  c(1, 2,
    3, 4,
    5, 6,
    7, 8),
  ncol = 2,
  byrow = TRUE
)

p_panel_8 <- grid.arrange(
  p_bor,
  p_rf,
  p_et,
  p_xgb,
  p_gbm,
  p_svm,
  p_heat,
  p_core_bar,
  layout_matrix = layout_mat
)

ggsave("Fig_AllModels_Heatmap_Core4_panel.pdf",
       p_panel_8, width = 50, height = 78, units = "cm", dpi = 400)

ROC：
library(readxl)
library(pROC)

f1 <- "C:/Users/Desktop/1.xlsx"
f2 <- "C:/Users/Desktop/2.xlsx"

metab <- read_excel(f1)
grp   <- read_excel(f2)
dat0  <- merge(metab, grp, by = "Sample")

orig_names <- colnames(dat0)

safe_names <- make.names(orig_names, unique = TRUE)
name_map   <- setNames(orig_names, safe_names)

colnames(dat0) <- safe_names
dat <- dat0

dat$Group_bin <- ifelse(dat$Group == "Sepsis", 1, 0)

feat_safe <- setdiff(colnames(dat), c("Sample", "Group", "Group_bin"))

roc_list <- lapply(feat_safe, function(x) roc(dat$Group_bin, dat[[x]]))
auc_list <- sapply(roc_list, function(r) as.numeric(auc(r)))

ord <- order(auc_list, decreasing = TRUE)
feat_safe <- feat_safe[ord]
roc_list  <- roc_list[ord]
auc_list  <- auc_list[ord]

cols <- grDevices::rainbow(length(feat_safe))

pdf("C:/Users/张天龙/Desktop/ROC_3x3_colored_sorted.pdf", width = 10, height = 10)
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))

for (i in seq_along(feat_safe)) {
  x_safe <- feat_safe[i]
  x_orig <- name_map[x_safe]
  
  r  <- roc_list[[i]]
  ci <- ci.auc(r, method = "delong")
  
  plot(r,
       col = cols[i], lwd = 2,
       legacy.axes = TRUE,
       main = paste0(
         x_orig,
         "\nAUC=", sprintf("%.2f", as.numeric(auc_list[i])),
         "  CI=", paste(sprintf("%.2f", as.numeric(ci)), collapse = "-")
       ))
  abline(a = 0, b = 1, lty = 2)
  grid()
}

dev.off()
cat("Generated: C:/Users/Desktop/ROC_3x3_colored_sorted.pdf\n")

Single CELL:
rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(harmony)
  library(clustree)
  library(SingleR)
  library(celldex)
  library(DoubletFinder)
  library(Matrix)
})

set.seed(42)
options(future.globals.maxSize = 6e9)

sample_paths <- c(
  GSM5333784 = "D:/data/single/analysis/GSM5333784/",
  GSM5333785 = "D:/data/single/analysis/GSM5333785/",
  GSM5333786 = "D:/data/single/analysis/GSM5333786/",
  GSM5333787 = "D:/data/single/analysis/GSM5333787/",
  GSM5333788 = "D:/data/single/analysis/GSM5333788/",
  GSM5333789 = "D:/data/single/analysis/GSM5333789/",
  GSM5333790 = "D:/data/single/analysis/GSM5333790/",
  GSM5333791 = "D:/data/single/analysis/GSM5333791/",
  GSM5333792 = "D:/data/single/analysis/GSM5333792/"
)

control_samples <- c(
  "GSM5333784", "GSM5333785",
  "GSM5333790", "GSM5333791", "GSM5333792"
)

out_base <- "D:/data/single/analysis/seurat_doubletfinder_final_results"
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_base, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_base, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_base, "DEG_results", "pseudobulk_DESeq2"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_base, "DEG_results", "single_cell_FindMarkers_fallback"), showWarnings = FALSE, recursive = TRUE)

missing_dirs <- names(sample_paths)[!dir.exists(sample_paths)]
if (length(missing_dirs) > 0) {
  stop("These sample folders do not exist: ", paste(missing_dirs, collapse = ", "))
}

mt_pattern <- "^MT-"
min_features <- 200
max_features <- 6000
max_percent_mt <- 10
pc.doublet <- 1:10
harmony.dims <- 1:30
cluster_resolution_for_annotation <- 0.3

get_df_fun <- function(...) {
  fn_names <- c(...)
  ns <- asNamespace("DoubletFinder")
  for (fn in fn_names) {
    if (exists(fn, envir = ns, mode = "function")) {
      return(get(fn, envir = ns))
    }
  }
  stop("Cannot find DoubletFinder function: ", paste(fn_names, collapse = " or "))
}

paramSweep_fun <- get_df_fun("paramSweep", "paramSweep_v3")
summarizeSweep_fun <- get_df_fun("summarizeSweep")
find_pK_fun <- get_df_fun("find.pK")
modelHomotypic_fun <- get_df_fun("modelHomotypic")
doubletFinder_fun <- get_df_fun("doubletFinder", "doubletFinder_v3")

read_10x_gene_expression <- function(data_dir) {
  x <- Read10X(data.dir = data_dir)
  if (is.list(x)) {
    if ("Gene Expression" %in% names(x)) {
      x <- x[["Gene Expression"]]
    } else {
      x <- x[[1]]
    }
  }
  x
}

get_assay_data_safe <- function(object, assay = "RNA", layer = "data") {
  tryCatch(
    GetAssayData(object, assay = assay, layer = layer),
    error = function(e) GetAssayData(object, assay = assay, slot = layer)
  )
}

to_broad_celltype <- function(x) {
  x <- ifelse(is.na(x), "", x)
  dplyr::case_when(
    grepl("platelet|megakaryo", x, ignore.case = TRUE) ~ "Platelet",
    grepl("eryth|red blood|rbc", x, ignore.case = TRUE) ~ "RBC",
    grepl("plasma", x, ignore.case = TRUE) ~ "B",
    grepl("natural killer|\\bnk\\b|nk_cell|nk cells", x, ignore.case = TRUE) ~ "NK",
    grepl("b[-_ ]?cell|b cells|\\bb_cell\\b", x, ignore.case = TRUE) ~ "B",
    grepl("t[-_ ]?cell|t cells|\\bt_cell\\b|cd4|cd8", x, ignore.case = TRUE) ~ "T",
    grepl("monocyte|macrophage", x, ignore.case = TRUE) ~ "Mono",
    grepl("dendritic|\\bdc\\b|\\bpdc\\b|\\bcdc\\b", x, ignore.case = TRUE) ~ "DC",
    grepl("neutrophil|granulocyte", x, ignore.case = TRUE) ~ "Neutrophil",
    grepl("progenitor|stem|hsc|cd34", x, ignore.case = TRUE) ~ "HSPC",
    TRUE ~ "Unknown"
  )
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_-]+", "_", x)
}

plot_volcano <- function(de, outfile, fc_col = "log2FoldChange", padj_col = "padj") {
  if (!all(c(fc_col, padj_col) %in% colnames(de))) return(invisible(NULL))
  de$neg_log10_padj <- -log10(de[[padj_col]] + 1e-300)
  de$status <- "NS"
  de$status[!is.na(de[[padj_col]]) & de[[padj_col]] < 0.05 & de[[fc_col]] > 0.25] <- "Up_in_Sep"
  de$status[!is.na(de[[padj_col]]) & de[[padj_col]] < 0.05 & de[[fc_col]] < -0.25] <- "Up_in_control"
  p <- ggplot(de, aes(x = .data[[fc_col]], y = neg_log10_padj, color = status)) +
    geom_point(size = 0.7, alpha = 0.75) +
    scale_color_manual(values = c(NS = "grey70", Up_in_Sep = "#D55E00", Up_in_control = "#0072B2")) +
    theme_classic() +
    labs(x = "log2FC: Sep vs control", y = "-log10 adjusted P", color = NULL)
  ggsave(outfile, p, width = 6.5, height = 5)
}

run_pseudobulk_deseq2 <- function(
    object,
    out_dir,
    celltype_col = "celltype",
    sample_col = "sample",
    group_col = "group",
    min_cells_per_sample = 10,
    min_samples_per_group = 2,
    min_count = 10,
    min_samples_with_count = 2
) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    message("DESeq2 is not installed; skip pseudobulk DESeq2.")
    return(NULL)
  }
  
  counts <- get_assay_data_safe(object, assay = "RNA", layer = "counts")
  meta <- object@meta.data
  celltypes <- sort(unique(as.character(meta[[celltype_col]])))
  celltypes <- setdiff(celltypes, c(NA, "", "Unknown", "RBC", "HSPC", "Cycling"))
  de_summary <- list()
  
  for (ct in celltypes) {
    message("Pseudobulk DESeq2: ", ct)
    cells <- rownames(meta)[as.character(meta[[celltype_col]]) == ct]
    sub_meta <- meta[cells, , drop = FALSE]
    
    sample_cell_counts <- sub_meta %>%
      dplyr::count(
        group = .data[[group_col]],
        sample = .data[[sample_col]],
        name = "n_cells"
      )
    
    write.csv(
      sample_cell_counts,
      file = file.path(out_dir, paste0(safe_name(ct), "_sample_coverage.csv")),
      row.names = FALSE
    )
    
    keep_samples <- sample_cell_counts$sample[sample_cell_counts$n_cells >= min_cells_per_sample]
    sub_meta <- sub_meta[sub_meta[[sample_col]] %in% keep_samples, , drop = FALSE]
    cells <- rownames(sub_meta)
    
    group_sample_counts <- sub_meta %>%
      dplyr::transmute(
        sample = .data[[sample_col]],
        group = .data[[group_col]]
      ) %>%
      dplyr::distinct() %>%
      dplyr::count(group, name = "n_samples")
    
    if (nrow(group_sample_counts) < 2 ||
        any(group_sample_counts$n_samples < min_samples_per_group)) {
      de_summary[[ct]] <- data.frame(celltype = ct, status = "skip_insufficient_samples")
      next
    }
    
    sample_factor <- factor(sub_meta[[sample_col]])
    design_mat <- Matrix::sparse.model.matrix(~ 0 + sample_factor)
    colnames(design_mat) <- levels(sample_factor)
    
    pb_counts <- counts[, cells, drop = FALSE] %*% design_mat
    pb_counts <- as.matrix(pb_counts)
    pb_counts <- round(pb_counts)
    storage.mode(pb_counts) <- "integer"
    
    coldata <- sub_meta %>%
      dplyr::transmute(
        sample = .data[[sample_col]],
        group = .data[[group_col]]
      ) %>%
      dplyr::distinct() %>%
      as.data.frame()
    rownames(coldata) <- coldata$sample
    coldata <- coldata[colnames(pb_counts), , drop = FALSE]
    coldata$group <- factor(coldata$group, levels = c("control", "Sep"))
    
    keep_genes <- rowSums(pb_counts >= min_count) >= min_samples_with_count
    pb_counts <- pb_counts[keep_genes, , drop = FALSE]
    
    if (nrow(pb_counts) < 100) {
      de_summary[[ct]] <- data.frame(celltype = ct, status = "skip_too_few_genes_after_filter")
      next
    }
    
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = pb_counts,
      colData = coldata,
      design = ~ group
    )
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds, contrast = c("group", "Sep", "control"))
    de <- as.data.frame(res)
    de$gene <- rownames(de)
    de <- de[order(de$padj, -abs(de$log2FoldChange), na.last = TRUE), ]
    
    write.csv(
      de,
      file = file.path(out_dir, paste0(safe_name(ct), "_pseudobulk_DESeq2_DEG.csv")),
      row.names = FALSE
    )
    
    plot_volcano(
      de,
      outfile = file.path(out_dir, paste0(safe_name(ct), "_pseudobulk_volcano.pdf")),
      fc_col = "log2FoldChange",
      padj_col = "padj"
    )
    
    de_summary[[ct]] <- data.frame(
      celltype = ct,
      status = "done",
      Sep_samples = sum(coldata$group == "Sep"),
      control_samples = sum(coldata$group == "control"),
      genes_tested = nrow(de),
      sig_genes_padj_0.05 = sum(!is.na(de$padj) & de$padj < 0.05)
    )
  }
  
  summary_df <- dplyr::bind_rows(de_summary)
  write.csv(summary_df, file.path(out_dir, "pseudobulk_DESeq2_summary.csv"), row.names = FALSE)
  summary_df
}

run_single_cell_findmarkers <- function(
    object,
    out_dir,
    celltype_col = "celltype",
    group_col = "group",
    min_cells_per_group = 20
) {
  celltypes <- sort(unique(as.character(object@meta.data[[celltype_col]])))
  celltypes <- setdiff(celltypes, c(NA, "", "Unknown", "RBC", "HSPC", "Cycling"))
  de_summary <- list()
  
  for (ct in celltypes) {
    message("Single-cell FindMarkers fallback: ", ct)
    ct_cells <- rownames(object@meta.data)[as.character(object@meta.data[[celltype_col]]) == ct]
    sub <- subset(object, cells = ct_cells)
    Idents(sub) <- group_col
    group_counts <- table(sub@meta.data[[group_col]])
    
    if (!all(c("Sep", "control") %in% names(group_counts)) ||
        any(group_counts[c("Sep", "control")] < min_cells_per_group)) {
      de_summary[[ct]] <- data.frame(celltype = ct, status = "skip_insufficient_cells")
      next
    }
    
    de <- FindMarkers(
      sub,
      ident.1 = "Sep",
      ident.2 = "control",
      min.pct = 0.1,
      logfc.threshold = 0.25
    )
    de$gene <- rownames(de)
    de <- de[order(de$p_val_adj, -abs(de$avg_log2FC), na.last = TRUE), ]
    
    write.csv(
      de,
      file = file.path(out_dir, paste0(safe_name(ct), "_FindMarkers_DEG.csv")),
      row.names = FALSE
    )
    
    plot_volcano(
      de,
      outfile = file.path(out_dir, paste0(safe_name(ct), "_FindMarkers_volcano.pdf")),
      fc_col = "avg_log2FC",
      padj_col = "p_val_adj"
    )
    
    de_summary[[ct]] <- data.frame(
      celltype = ct,
      status = "done",
      Sep_cells = as.numeric(group_counts["Sep"]),
      control_cells = as.numeric(group_counts["control"]),
      genes_tested = nrow(de),
      sig_genes_padj_0.05 = sum(!is.na(de$p_val_adj) & de$p_val_adj < 0.05)
    )
  }
  
  summary_df <- dplyr::bind_rows(de_summary)
  write.csv(summary_df, file.path(out_dir, "FindMarkers_summary.csv"), row.names = FALSE)
  summary_df
}

old_assay_opt <- getOption("Seurat.object.assay.version")
options(Seurat.object.assay.version = "v3")

singlet_mats <- list()
qc_summary <- list()
doublet_summary <- list()
pk_summary <- list()

for (sid in names(sample_paths)) {
  message("Processing sample: ", sid)
  
  mat <- read_10x_gene_expression(sample_paths[[sid]])
  raw_barcodes <- ncol(mat)
  
  seu_tmp <- CreateSeuratObject(
    counts = mat,
    project = sid,
    min.cells = 3,
    min.features = min_features
  )
  cells_after_initial_filter <- ncol(seu_tmp)
  
  seu_tmp$sample <- sid
  seu_tmp$group <- ifelse(sid %in% control_samples, "control", "Sep")
  seu_tmp[["percent.mt"]] <- PercentageFeatureSet(seu_tmp, pattern = mt_pattern)
  
  seu_tmp <- subset(
    seu_tmp,
    subset = nFeature_RNA > min_features &
      nFeature_RNA < max_features &
      percent.mt < max_percent_mt
  )
  qc_after <- ncol(seu_tmp)
  if (qc_after < 100) stop("Too few cells after QC for ", sid, ": ", qc_after)
  
  qc_summary[[sid]] <- data.frame(
    sample = sid,
    raw_barcodes = raw_barcodes,
    cells_after_initial_min_features = cells_after_initial_filter,
    cells_after_qc = qc_after,
    group = unique(seu_tmp$group)
  )
  
  seu_tmp <- NormalizeData(seu_tmp, verbose = FALSE)
  seu_tmp <- FindVariableFeatures(seu_tmp, nfeatures = 2000, verbose = FALSE)
  seu_tmp <- ScaleData(seu_tmp, verbose = FALSE)
  seu_tmp <- RunPCA(seu_tmp, npcs = max(pc.doublet), seed.use = 42, verbose = FALSE)
  seu_tmp <- RunUMAP(seu_tmp, dims = pc.doublet, seed.use = 42, verbose = FALSE)
  seu_tmp <- FindNeighbors(seu_tmp, dims = pc.doublet, verbose = FALSE)
  seu_tmp <- FindClusters(seu_tmp, resolution = 0.3, random.seed = 42, verbose = FALSE)
  
  sweep.res <- paramSweep_fun(seu_tmp, PCs = pc.doublet, sct = FALSE)
  sweep.stats <- summarizeSweep_fun(sweep.res, GT = FALSE)
  bcmvn <- find_pK_fun(sweep.stats)
  bcmvn <- bcmvn[order(-bcmvn$BCmetric), , drop = FALSE]
  opt_pK <- as.numeric(as.character(bcmvn$pK[1]))
  message("Best pK for ", sid, " = ", opt_pK)
  
  homotypic.prop <- tryCatch(
    modelHomotypic_fun(as.character(Idents(seu_tmp))),
    error = function(e) {
      message("modelHomotypic failed for ", sid, "; use no homotypic correction.")
      0
    }
  )
  
  doublet_rate_sid <- min(0.10, ncol(seu_tmp) * 8e-6)
  nExp_poi <- round(doublet_rate_sid * ncol(seu_tmp))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  nExp_poi.adj <- max(1, min(nExp_poi.adj, ncol(seu_tmp) - 1))
  
  pk_summary[[sid]] <- data.frame(
    sample = sid,
    pK = opt_pK,
    BCmetric = bcmvn$BCmetric[1],
    doublet_rate = doublet_rate_sid,
    homotypic_prop = homotypic.prop,
    nExp_raw = nExp_poi,
    nExp_homotypic_adjusted = nExp_poi.adj
  )
  
  seu_tmp <- doubletFinder_fun(
    seu_tmp,
    PCs = pc.doublet,
    pN = 0.25,
    pK = opt_pK,
    nExp = nExp_poi.adj,
    reuse.pANN = NULL,
    sct = FALSE
  )
  
  df_col <- grep("^DF.classifications", colnames(seu_tmp@meta.data), value = TRUE)
  df_col <- tail(df_col, 1)
  if (length(df_col) == 0) stop("Cannot find DoubletFinder classification column for ", sid)
  
  df_tab <- table(seu_tmp@meta.data[[df_col]])
  doublet_summary[[sid]] <- data.frame(
    sample = sid,
    classification = names(df_tab),
    n_cells = as.numeric(df_tab)
  )
  
  pdf(file.path(out_base, "plots", paste0(sid, "_DoubletFinder_UMAP.pdf")), width = 7, height = 5)
  print(DimPlot(seu_tmp, group.by = df_col))
  dev.off()
  
  pdf(file.path(out_base, "plots", paste0(sid, "_QC_by_DoubletFinder.pdf")), width = 10, height = 4)
  print(VlnPlot(seu_tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = df_col, ncol = 3, pt.size = 0))
  dev.off()
  
  singlet_cells <- rownames(seu_tmp@meta.data)[seu_tmp@meta.data[[df_col]] == "Singlet"]
  singlet_mats[[sid]] <- mat[, singlet_cells, drop = FALSE]
}

write.csv(dplyr::bind_rows(qc_summary), file.path(out_base, "tables", "QC_summary.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(pk_summary), file.path(out_base, "tables", "DoubletFinder_best_pK.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(doublet_summary), file.path(out_base, "tables", "DoubletFinder_summary.csv"), row.names = FALSE)

if (is.null(old_assay_opt)) {
  options(Seurat.object.assay.version = "v5")
} else {
  options(Seurat.object.assay.version = old_assay_opt)
}

obj_list <- lapply(names(singlet_mats), function(sid) {
  seu <- CreateSeuratObject(counts = singlet_mats[[sid]], project = sid, min.cells = 3)
  seu$sample <- sid
  seu$group <- ifelse(sid %in% control_samples, "control", "Sep")
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = mt_pattern)
  seu
})
names(obj_list) <- names(singlet_mats)

merged <- merge(
  x = obj_list[[1]],
  y = obj_list[-1],
  add.cell.ids = names(obj_list),
  project = "PBMC_9samples"
)

merged$sample <- merged$orig.ident
merged$group <- ifelse(merged$sample %in% control_samples, "control", "Sep")
merged$group <- factor(merged$group, levels = c("control", "Sep"))

write.csv(
  dplyr::count(merged[[]], sample, group, name = "n_cells_after_doublet_filter"),
  file.path(out_base, "tables", "cell_counts_after_doublet_filter.csv"),
  row.names = FALSE
)

merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, nfeatures = 2000, verbose = FALSE)
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, features = VariableFeatures(merged), npcs = max(harmony.dims), seed.use = 42, verbose = FALSE)

pdf(file.path(out_base, "plots", "ElbowPlot.pdf"), width = 7, height = 5)
print(ElbowPlot(merged, ndims = 50))
dev.off()

merged <- RunHarmony(
  object = merged,
  group.by.vars = "sample",
  reduction = "pca",
  dims.use = harmony.dims,
  reduction.save = "harmony",
  verbose = FALSE
)

merged <- FindNeighbors(merged, reduction = "harmony", dims = harmony.dims, verbose = FALSE)

resolution_values <- seq(0.1, 0.9, by = 0.1)
for (res in resolution_values) {
  merged <- FindClusters(
    merged,
    resolution = res,
    cluster.name = paste0("harmony_clusters_", res),
    random.seed = 42,
    verbose = FALSE
  )
}

pdf(file.path(out_base, "plots", "clustree_harmony.pdf"), width = 12, height = 10)
print(clustree(merged@meta.data, prefix = "harmony_clusters_"))
dev.off()

cluster_col <- paste0("harmony_clusters_", cluster_resolution_for_annotation)
if (!cluster_col %in% colnames(merged@meta.data)) {
  stop("Cannot find cluster column: ", cluster_col)
}

Idents(merged) <- cluster_col
merged$harmony_clusters <- as.character(Idents(merged))

merged <- RunUMAP(merged, reduction = "harmony", dims = harmony.dims, seed.use = 42, verbose = FALSE)

pdf(file.path(out_base, "plots", "UMAP_harmony_clusters.pdf"), width = 8, height = 6)
print(DimPlot(merged, group.by = "harmony_clusters", label = TRUE))
dev.off()

pdf(file.path(out_base, "plots", "UMAP_group.pdf"), width = 8, height = 6)
print(DimPlot(merged, group.by = "group"))
dev.off()

pdf(file.path(out_base, "plots", "UMAP_sample.pdf"), width = 9, height = 6)
print(DimPlot(merged, group.by = "sample"))
dev.off()

if (inherits(merged[["RNA"]], "Assay5") && length(Layers(merged[["RNA"]])) > 1) {
  merged <- JoinLayers(merged)
}
DefaultAssay(merged) <- "RNA"

refdata <- tryCatch(
  celldex::MonacoImmuneData(),
  error = function(e) {
    message("MonacoImmuneData failed; use HumanPrimaryCellAtlasData instead.")
    celldex::HumanPrimaryCellAtlasData()
  }
)

ref_coldata <- SummarizedExperiment::colData(refdata)
main_labels <- if ("label.main" %in% colnames(ref_coldata)) refdata$label.main else refdata$label.fine
fine_labels <- if ("label.fine" %in% colnames(ref_coldata)) refdata$label.fine else refdata$label.main

test_mat <- get_assay_data_safe(merged, assay = "RNA", layer = "data")

pred_main <- SingleR(
  test = test_mat,
  ref = refdata,
  labels = main_labels,
  clusters = merged$harmony_clusters,
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts"
)

pred_fine <- SingleR(
  test = test_mat,
  ref = refdata,
  labels = fine_labels,
  clusters = merged$harmony_clusters,
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts"
)

pdf(file.path(out_base, "plots", "SingleR_main_heatmap.pdf"), width = 12, height = 7)
plotScoreHeatmap(pred_main, cellwidth = 20, cellheight = 9)
dev.off()

pdf(file.path(out_base, "plots", "SingleR_fine_heatmap.pdf"), width = 14, height = 8)
plotScoreHeatmap(pred_fine, cellwidth = 20, cellheight = 9)
dev.off()

pred_df <- data.frame(
  harmony_clusters = rownames(pred_fine),
  SingleR_main = pred_main$labels[match(rownames(pred_fine), rownames(pred_main))],
  SingleR_fine = pred_fine$labels,
  SingleR_pruned = pred_fine$pruned.labels,
  stringsAsFactors = FALSE
)
pred_df$SingleR_label <- ifelse(is.na(pred_df$SingleR_pruned), pred_df$SingleR_fine, pred_df$SingleR_pruned)
pred_df$SingleR_broad <- to_broad_celltype(pred_df$SingleR_label)

merged$SingleR_main <- pred_df$SingleR_main[match(merged$harmony_clusters, pred_df$harmony_clusters)]
merged$SingleR_label <- pred_df$SingleR_label[match(merged$harmony_clusters, pred_df$harmony_clusters)]

write.csv(pred_df, file.path(out_base, "tables", "SingleR_cluster_labels.csv"), row.names = FALSE)

canonical_markers <- list(
  T = c("CD3D", "CD3E", "TRAC", "IL7R", "CCR7", "LTB", "CD4", "CD8A", "CD8B"),
  NK = c("NKG7", "GNLY", "PRF1", "GZMB", "KLRD1", "KLRF1", "FCGR3A"),
  B = c("MS4A1", "CD79A", "CD79B", "CD74", "CD37", "BANK1"),
  Plasma = c("MZB1", "JCHAIN", "IGHG1", "IGKC"),
  Mono = c("LYZ", "LST1", "S100A8", "S100A9", "FCN1", "CD14", "VCAN", "FCGR3A", "MS4A7"),
  DC = c("FCER1A", "CLEC10A", "CD1C", "CST3", "LILRA4", "CLEC4C", "TCF4", "IRF7"),
  Neutrophil = c("CSF3R", "S100A12", "LTF", "MMP9", "FCGR3B", "CEACAM8", "CXCR2"),
  Platelet = c("PPBP", "PF4", "GP9", "NRGN"),
  RBC = c("HBB", "HBA1", "HBA2", "ALAS2", "AHSP", "GYPB"),
  Cycling = c("MKI67", "TOP2A", "STMN1", "TYMS", "RRM2"),
  HSPC = c("CD34", "PROM1", "SPINK2", "AVP")
)

markers.to.plot <- intersect(unique(unlist(canonical_markers)), rownames(merged))

pdf(file.path(out_base, "plots", "DotPlot_canonical_markers.pdf"), width = 22, height = 9)
print(DotPlot(merged, features = markers.to.plot, group.by = "harmony_clusters") + RotatedAxis())
dev.off()

score_cols <- c()
score_names <- c()
for (ct in names(canonical_markers)) {
  feats <- intersect(canonical_markers[[ct]], rownames(merged))
  if (length(feats) >= 2) {
    old_cols <- colnames(merged[[]])
    merged <- AddModuleScore(
      object = merged,
      features = list(feats),
      name = paste0("score_", ct, "_"),
      assay = "RNA"
    )
    new_cols <- setdiff(colnames(merged[[]]), old_cols)
    score_cols <- c(score_cols, new_cols[1])
    score_names <- c(score_names, ct)
  }
}

cluster_score_df <- merged[[]] %>%
  dplyr::select(harmony_clusters, dplyr::all_of(score_cols)) %>%
  dplyr::group_by(harmony_clusters) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(score_cols), mean), .groups = "drop")

score_only <- as.data.frame(cluster_score_df[, score_cols, drop = FALSE])
colnames(score_only) <- score_names
cluster_score_df$marker_top <- colnames(score_only)[max.col(score_only, ties.method = "first")]

write.csv(cluster_score_df, file.path(out_base, "tables", "cluster_marker_module_scores.csv"), row.names = FALSE)

Idents(merged) <- "harmony_clusters"
cluster_markers <- FindAllMarkers(
  merged,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

write.csv(cluster_markers, file.path(out_base, "tables", "cluster_markers_all.csv"), row.names = FALSE)

top15_markers <- cluster_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 15, with_ties = FALSE) %>%
  dplyr::ungroup()

write.csv(top15_markers, file.path(out_base, "tables", "cluster_markers_top15.csv"), row.names = FALSE)

annotation_evidence <- top15_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(top_genes = paste(gene, collapse = ", "), .groups = "drop") %>%
  dplyr::mutate(harmony_clusters = as.character(cluster)) %>%
  dplyr::select(harmony_clusters, top_genes)

annotation_table <- pred_df %>%
  dplyr::left_join(cluster_score_df %>% dplyr::select(harmony_clusters, marker_top), by = "harmony_clusters") %>%
  dplyr::left_join(annotation_evidence, by = "harmony_clusters") %>%
  dplyr::arrange(as.numeric(harmony_clusters))

suggested <- annotation_table$marker_top
suggested[suggested == "Plasma"] <- "B"
bad_suggested <- is.na(suggested) | suggested == ""
suggested[bad_suggested] <- annotation_table$SingleR_broad[bad_suggested]
suggested[suggested %in% c("RBC", "HSPC", "Cycling", "Unknown")] <- "Unknown"
annotation_table$celltype_suggested <- suggested
annotation_table$needs_review <- annotation_table$SingleR_broad != annotation_table$marker_top

write.csv(annotation_table, file.path(out_base, "tables", "cluster_annotation_evidence.csv"), row.names = FALSE)

manual_celltype_map <- c(
  "0"  = "Mono",
  "1"  = "Mono",
  "2"  = "Neutrophil",
  "3"  = "T",
  "4"  = "Platelet",
  "5"  = "NK",
  "6"  = "Mono",
  "7"  = "T",
  "8"  = "RBC",
  "9"  = "DC",
  "10" = "T",
  "11" = "B",
  "12" = "Plasma",
  "13" = "DC",
  "14" = "Neutrophil",
  "15" = "NK",
  "16" = "Unknown"
)

cluster_to_celltype <- manual_celltype_map

missing_clusters <- setdiff(unique(as.character(merged$harmony_clusters)), names(cluster_to_celltype))
if (length(missing_clusters) > 0) {
  stop("manual_celltype_map is missing these clusters: ", paste(missing_clusters, collapse = ", "))
}

merged$celltype <- unname(cluster_to_celltype[as.character(merged$harmony_clusters)])
merged$celltype[is.na(merged$celltype)] <- "Unknown"

merged$celltype <- factor(
  merged$celltype,
  levels = c(
    "Mono", "T", "NK", "B", "Plasma",
    "DC", "Neutrophil", "Platelet",
    "RBC", "Unknown"
  )
)

write.csv(
  dplyr::count(merged[[]], harmony_clusters, celltype, group, sample, name = "n_cells"),
  file.path(out_base, "tables", "celltype_group_sample_counts_all.csv"),
  row.names = FALSE
)

pdf(file.path(out_base, "plots", "UMAP_celltype_manual.pdf"), width = 10, height = 7)
print(DimPlot(merged, group.by = "celltype", label = TRUE, repel = TRUE))
dev.off()

pdf(file.path(out_base, "plots", "DotPlot_canonical_markers_by_manual_celltype.pdf"), width = 22, height = 9)
print(DotPlot(merged, features = markers.to.plot, group.by = "celltype") + RotatedAxis())
dev.off()

deg_object <- subset(
  merged,
  subset = celltype %in% c(
    "Mono", "T", "NK", "B", "Plasma",
    "DC", "Neutrophil", "Platelet"
  )
)

Idents(deg_object) <- "celltype"

write.csv(
  dplyr::count(deg_object[[]], celltype, group, sample, name = "n_cells"),
  file.path(out_base, "tables", "celltype_group_sample_counts_for_DEG.csv"),
  row.names = FALSE
)

pseudobulk_summary <- run_pseudobulk_deseq2(
  object = deg_object,
  out_dir = file.path(out_base, "DEG_results", "pseudobulk_DESeq2")
)

if (is.null(pseudobulk_summary)) {
  run_single_cell_findmarkers(
    object = deg_object,
    out_dir = file.path(out_base, "DEG_results", "single_cell_FindMarkers_fallback")
  )
}

dir.create(
  file.path(out_base, "DEG_results", "single_cell_FindMarkers_fallback"),
  showWarnings = FALSE,
  recursive = TRUE
)

run_single_cell_findmarkers(
  object = deg_object,
  out_dir = file.path(out_base, "DEG_results", "single_cell_FindMarkers_fallback")
)

saveRDS(
  merged,
  file = file.path(out_base, "merged_harmony_doubletfinder_annotated.rds")
)

saveRDS(
  deg_object,
  file = file.path(out_base, "merged_for_celltype_DEG.rds")
)

message("All done. Results directory: ", out_base)
message("Please check annotation evidence first: ", file.path(out_base, "tables", "cluster_annotation_evidence.csv"))
