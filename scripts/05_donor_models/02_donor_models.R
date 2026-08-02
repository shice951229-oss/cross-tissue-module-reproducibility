########################################
# Script: 02_donor_models.R
# Purpose: Recover verified GSE102556 donors and compare donor-aware and naive models for GSE102556 and GSE48350.
# Input: Locked module scores, phenotype metadata, and verified donor maps.
# Output: Donor crosswalks, model result tables, and comparison figures.
# Software: R
# Version: 4.5.2
# Random seed: Not applicable (deterministic)
# Author: Study authors
########################################
options(stringsAsFactors = FALSE, width = 180, contrasts = c("contr.treatment", "contr.poly"))

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(ggplot2)
  library(patchwork)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out <- file.path(root, "13_methodological_reanalysis_v9")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(20260802)

write_csv <- function(x, name) utils::write.csv(x, file.path(out, name), row.names = FALSE, na = "")
num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.+-]", "", as.character(x))))
z <- function(x) as.numeric(scale(x))
fdr_status <- function(x) ifelse(is.finite(x) & x < 0.05, "FDR<0.05", "not_FDR<0.05")

extract_fixed <- function(fit, term) {
  s <- coef(summary(fit))
  if (!term %in% rownames(s)) return(c(beta = NA, SE = NA, df = NA, CI_low = NA, CI_high = NA, P = NA))
  pcol <- grep("Pr\\(", colnames(s), value = TRUE)[1]
  df <- if ("df" %in% colnames(s)) s[term, "df"] else Inf
  crit <- if (is.finite(df)) qt(.975, df) else qnorm(.975)
  c(beta = s[term, "Estimate"], SE = s[term, "Std. Error"], df = df,
    CI_low = s[term, "Estimate"] - crit * s[term, "Std. Error"],
    CI_high = s[term, "Estimate"] + crit * s[term, "Std. Error"],
    P = s[term, pcol])
}

rank_ok <- function(formula, data) {
  mm <- model.matrix(formula, data = data)
  qr(mm)$rank == ncol(mm)
}

save_plot <- function(p, stem, width_mm = 183, height_mm = 120) {
  ggplot2::ggsave(file.path(out, paste0(stem, ".pdf")), p, width = width_mm / 25.4,
                  height = height_mm / 25.4, device = grDevices::cairo_pdf, family = "Arial")
  ggplot2::ggsave(file.path(out, paste0(stem, ".png")), p, width = width_mm / 25.4,
                  height = height_mm / 25.4, dpi = 600, device = ragg::agg_png)
}

theme_set(theme_classic(base_size = 7, base_family = "Arial") +
            theme(axis.line = element_line(linewidth = .35), axis.ticks = element_line(linewidth = .35),
                  legend.title = element_text(size = 6.5), legend.text = element_text(size = 6),
                  strip.text = element_text(size = 6.5, face = "bold"),
                  plot.title = element_text(size = 8, face = "bold"),
                  plot.tag = element_text(size = 9, face = "bold")))

# -------------------------------------------------------------------------
# GSE102556 donor crosswalk and donor-aware models
# -------------------------------------------------------------------------
ph <- readRDS(file.path(root, "04_intermediate", "phenotype_tables", "MDD_GSE102556", "pheno_harmonized.rds"))
scores <- read.csv(file.path(root, "05_results", "MDD_projection", "MDD_brain_projection_module_scores.csv"), check.names = FALSE)
module_names <- names(scores)[2:25]

title_donor <- trimws(sub(":.*$", "", as.character(ph$title)))
sample_donor <- sub("\\..*$", "", as.character(ph$sample_id))
fpkm_header <- names(utils::read.delim(gzfile(file.path(root, "01_raw_data", "depression_MDD", "supplementary", "GSE102556", "GSE102556_HumanMDD_fpkmtab.txt.gz")), nrows = 1L, check.names = FALSE))
fpkm_samples <- setdiff(fpkm_header, c("gene_id", "gene_name", "biotype", "length"))

within_conflict <- function(v, donor) {
  ave(as.character(v), donor, FUN = function(y) length(unique(y[!is.na(y) & nzchar(y)])) > 1L) == "TRUE"
}

cw <- data.frame(
  GSM = as.character(ph$geo_accession),
  donor_id = as.character(ph$subject_id),
  diagnosis = as.character(ph$diagnosis),
  sex = as.character(ph[["gender:ch1"]]),
  brain_region = as.character(ph$brain_region),
  age = num(ph[["age:ch1"]]),
  RIN = num(ph[["rin:ch1"]]),
  PMI = num(ph[["pmi:ch1"]]),
  batch = NA_character_,
  instrument_model = as.character(ph$instrument_model),
  library_strategy = as.character(ph$library_strategy),
  sample_title = as.character(ph$title),
  source_name = as.character(ph$source_name_ch1),
  sample_id = as.character(ph$sample_id),
  donor_id_source = "GEO sample title explicit prefix before ':'; corroborated by project sample_id and FPKM column name",
  title_donor_id = title_donor,
  sample_id_donor_id = sample_donor,
  fpkm_column_present = as.character(ph$sample_id) %in% fpkm_samples,
  stringsAsFactors = FALSE
)
cw$conflict_donor_fields <- cw$donor_id != cw$title_donor_id | cw$donor_id != cw$sample_id_donor_id
cw$conflict_within_donor_diagnosis <- within_conflict(cw$diagnosis, cw$donor_id)
cw$conflict_within_donor_sex <- within_conflict(cw$sex, cw$donor_id)
cw$conflict_within_donor_age <- within_conflict(cw$age, cw$donor_id)
cw$conflict <- cw$conflict_donor_fields | !cw$fpkm_column_present |
  cw$conflict_within_donor_diagnosis | cw$conflict_within_donor_sex | cw$conflict_within_donor_age
cw$donor_id_verification_status <- ifelse(is.na(cw$donor_id) | !nzchar(cw$donor_id), "missing",
                                           ifelse(cw$conflict, "conflict", "verified_by_three_project_fields"))
write_csv(cw, "GSE102556_sample_donor_crosswalk.csv")

md <- merge(scores, cw, by.x = "sample_id", by.y = "sample_id", all.x = TRUE, sort = FALSE,
            suffixes = c("_score", ""))
if (nrow(md) != nrow(scores) || anyNA(md$donor_id)) stop("GSE102556 score/crosswalk alignment failed")
md$diagnosis <- factor(md$diagnosis, levels = c("Control", "MDD"))
md$sex <- factor(tolower(md$sex), levels = c("male", "female"))
md$brain_region <- factor(md$brain_region)
md$donor_id <- factor(md$donor_id)
for (nm in c("age", "RIN", "PMI")) md[[paste0(nm, "_z")]] <- z(md[[nm]])

cov_audit <- do.call(rbind, lapply(c("age", "sex", "RIN", "PMI"), function(v) {
  value <- md[[v]]
  data.frame(covariate = v, missing_n = sum(is.na(value)), missing_fraction = mean(is.na(value)),
             unique_nonmissing = length(unique(value[!is.na(value)])),
             values_by_diagnosis = paste(capture.output(print(tapply(value, md$diagnosis, function(x) {
               if (is.numeric(x)) sprintf("n=%d mean=%.3f sd=%.3f", sum(is.finite(x)), mean(x, na.rm=TRUE), sd(x, na.rm=TRUE)) else paste(table(x), collapse=";")
             }))), collapse = " "), stringsAsFactors = FALSE)
}))

base_vars_mdd <- c("diagnosis", "brain_region", "donor_id", "age_z", "sex", "RIN_z", "PMI_z")
analytic_mdd <- md[complete.cases(md[, base_vars_mdd]), ]
fixed_formula_mdd <- score_z ~ diagnosis + brain_region + age_z + sex + RIN_z + PMI_z
if (!rank_ok(fixed_formula_mdd, transform(analytic_mdd, score_z = 0))) stop("GSE102556 fixed-effects design is rank deficient")

fit_mdd_module <- function(module) {
  d <- analytic_mdd
  d$score_z <- z(d[[module]])
  fit_da <- lmer(score_z ~ diagnosis + brain_region + age_z + sex + RIN_z + PMI_z + (1 | donor_id),
                 data = d, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
  fit_nv <- lm(score_z ~ diagnosis + brain_region + age_z + sex + RIN_z + PMI_z, data = d)
  e_da <- extract_fixed(fit_da, "diagnosisMDD")
  e_nv <- extract_fixed(fit_nv, "diagnosisMDD")
  vc <- as.data.frame(VarCorr(fit_da))
  data.frame(module = module,
             donoraware_beta = e_da["beta"], donoraware_SE = e_da["SE"], donoraware_df = e_da["df"],
             donoraware_CI_low = e_da["CI_low"], donoraware_CI_high = e_da["CI_high"], donoraware_P = e_da["P"],
             naive_beta = e_nv["beta"], naive_SE = e_nv["SE"], naive_df = e_nv["df"],
             naive_CI_low = e_nv["CI_low"], naive_CI_high = e_nv["CI_high"], naive_P = e_nv["P"],
             sample_n = nrow(d), donor_n = nlevels(d$donor_id),
             random_intercept_variance = vc$vcov[vc$grp == "donor_id"][1],
             residual_variance = vc$vcov[vc$grp == "Residual"][1],
             ICC = vc$vcov[vc$grp == "donor_id"][1] / sum(vc$vcov[vc$grp %in% c("donor_id", "Residual")]),
             singular_fit = isSingular(fit_da, tol = 1e-5),
             donoraware_formula = "score_z ~ diagnosis + brain_region + age_z + sex + RIN_z + PMI_z + (1 | donor_id)",
             naive_formula = "score_z ~ diagnosis + brain_region + age_z + sex + RIN_z + PMI_z",
             stringsAsFactors = FALSE)
}
mdd_res <- do.call(rbind, lapply(module_names, fit_mdd_module))
mdd_res$donoraware_BH_FDR <- p.adjust(mdd_res$donoraware_P, "BH")
mdd_res$naive_BH_FDR <- p.adjust(mdd_res$naive_P, "BH")
mdd_res$beta_change <- mdd_res$donoraware_beta - mdd_res$naive_beta
mdd_res$SE_change <- mdd_res$donoraware_SE - mdd_res$naive_SE
mdd_res$naive_SE_over_donoraware_SE <- mdd_res$naive_SE / mdd_res$donoraware_SE
mdd_res$significance_transition <- paste(fdr_status(mdd_res$naive_BH_FDR), "to", fdr_status(mdd_res$donoraware_BH_FDR))
write_csv(mdd_res, "GSE102556_naive_vs_donoraware.csv")
write_csv(mdd_res[, c("module", grep("^donoraware_|^sample_n$|^donor_n$|variance$|ICC$|singular_fit$", names(mdd_res), value=TRUE))], "GSE102556_overall_model_results.csv")

# Sex interaction, pre-specified sensitivity because the source study was sex-stratified.
interaction_rows <- lapply(module_names, function(module) {
  d <- analytic_mdd; d$score_z <- z(d[[module]])
  fit <- lmer(score_z ~ diagnosis * sex + brain_region + age_z + RIN_z + PMI_z + (1 | donor_id),
              data = d, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
  term <- grep("diagnosisMDD:sex", rownames(coef(summary(fit))), value = TRUE)[1]
  e <- extract_fixed(fit, term)
  data.frame(module = module, interaction_term = term, beta = e["beta"], SE = e["SE"], df = e["df"],
             CI_low = e["CI_low"], CI_high = e["CI_high"], P = e["P"], sample_n=nrow(d), donor_n=nlevels(d$donor_id),
             formula = "score_z ~ diagnosis * sex + brain_region + age_z + RIN_z + PMI_z + (1 | donor_id)")
})
interaction_res <- do.call(rbind, interaction_rows)
interaction_res$BH_FDR <- p.adjust(interaction_res$P, "BH")
write_csv(interaction_res, "GSE102556_sex_interaction_sensitivity.csv")

# Region-specific descriptive family: global 144-test BH plus within-region 24-test BH.
region_rows <- list()
for (reg in levels(md$brain_region)) {
  dr <- md[md$brain_region == reg, ]
  duplicate_in_region <- any(table(dr$donor_id) > 1L)
  for (module in module_names) {
    vars <- c("diagnosis", "age_z", "sex", "RIN_z", "PMI_z", if (duplicate_in_region) "donor_id" else character())
    d <- droplevels(dr[complete.cases(dr[, vars]), ])
    d$score_z <- z(d[[module]])
    fixed_terms <- c("diagnosis", "age_z", "sex", "RIN_z", "PMI_z")
    candidate <- as.formula(paste("score_z ~", paste(fixed_terms, collapse = " + ")))
    dropped <- character()
    for (drop_nm in c("PMI_z", "RIN_z", "age_z", "sex")) {
      if (rank_ok(candidate, d)) break
      fixed_terms <- setdiff(fixed_terms, drop_nm); dropped <- c(dropped, drop_nm)
      candidate <- as.formula(paste("score_z ~", paste(fixed_terms, collapse = " + ")))
    }
    if (!rank_ok(candidate, d)) stop("Unresolved rank deficiency in GSE102556 region model: ", reg)
    if (duplicate_in_region) {
      full_formula <- as.formula(paste(deparse(candidate), "+ (1 | donor_id)"))
      fit <- lmer(full_formula, data=d, REML=TRUE, control=lmerControl(optimizer="bobyqa"))
      singular <- isSingular(fit, tol=1e-5)
    } else {
      full_formula <- candidate; fit <- lm(full_formula, data=d); singular <- FALSE
    }
    e <- extract_fixed(fit, "diagnosisMDD")
    region_rows[[length(region_rows)+1L]] <- data.frame(
      brain_region=reg, module=module, beta=e["beta"], SE=e["SE"], df=e["df"], CI_low=e["CI_low"], CI_high=e["CI_high"], P=e["P"],
      sample_n=nrow(d), donor_n=length(unique(d$donor_id)), MDD_sample_n=sum(d$diagnosis=="MDD"), control_sample_n=sum(d$diagnosis=="Control"),
      duplicate_donor_within_region=duplicate_in_region, singular_fit=singular,
      formula=paste(deparse(full_formula), collapse=" "), dropped_covariates=paste(dropped,collapse=";"), stringsAsFactors=FALSE)
  }
}
region_res <- do.call(rbind, region_rows)
region_res$BH_FDR_global_144 <- p.adjust(region_res$P, "BH")
region_res$BH_FDR_within_region_24 <- ave(region_res$P, region_res$brain_region, FUN=function(x)p.adjust(x,"BH"))
write_csv(region_res, "GSE102556_region_results.csv")

# Donor audit report
donor_counts <- aggregate(brain_region ~ donor_id, cw, function(x) length(unique(x)))
names(donor_counts)[2] <- "n_regions"
region_summary <- aggregate(GSM ~ brain_region, cw, length); names(region_summary)[2] <- "samples"
region_donors <- aggregate(donor_id ~ brain_region, cw, function(x) length(unique(x))); names(region_donors)[2] <- "donors"
region_summary <- merge(region_summary, region_donors, by="brain_region")
dx_donors <- aggregate(donor_id ~ diagnosis, cw, function(x) length(unique(x))); names(dx_donors)[2] <- "donors"
report <- c(
  "# GSE102556 donor audit report", "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")), "",
  "## Identity evidence", "",
  "The donor identifier is the explicit numeric prefix before the colon in each GEO sample title. It was cross-checked against the project `subject_id`, the prefix of the expression-table sample name, and presence of the complete sample name among the downloaded FPKM columns. No donor was inferred from ordering or diagnosis.", "",
  paste0("- Total samples: ", nrow(cw)),
  paste0("- Independent donors: ", length(unique(cw$donor_id))),
  paste0("- Missing donor identifiers: ", sum(is.na(cw$donor_id) | !nzchar(cw$donor_id))),
  paste0("- Samples with any donor conflict: ", sum(cw$conflict)),
  paste0("- Donor-aware analytic samples: ", nrow(analytic_mdd)), "",
  "## Donors by diagnosis", "",
  paste(capture.output(print(dx_donors, row.names=FALSE)), collapse="\n"), "",
  "## Regions", "",
  paste(capture.output(print(region_summary, row.names=FALSE)), collapse="\n"), "",
  "## Regions contributed per donor", "",
  paste(capture.output(print(table(donor_counts$n_regions))), collapse="\n"), "",
  "## Covariate audit", "",
  paste(capture.output(print(cov_audit, row.names=FALSE)), collapse="\n"), "",
  "## Primary model", "",
  "`score_z ~ diagnosis + brain_region + age_z + sex + RIN_z + PMI_z + (1 | donor_id)`", "",
  "The 24 overall diagnosis effects form one BH family. The naive comparison uses the identical analytic sample set and fixed-effect design but omits the donor random intercept. Region-specific analyses are descriptive; their primary correction spans all 144 tests, with within-region 24-test BH values supplied only as sensitivity results.", "",
  paste0("Effect correlation (naive vs donor-aware): ", sprintf("%.4f", cor(mdd_res$naive_beta, mdd_res$donoraware_beta, use="complete.obs"))),
  paste0("Naive FDR<0.05 modules: ", sum(mdd_res$naive_BH_FDR < .05, na.rm=TRUE)),
  paste0("Donor-aware FDR<0.05 modules: ", sum(mdd_res$donoraware_BH_FDR < .05, na.rm=TRUE)),
  paste0("Median naive SE / donor-aware SE: ", sprintf("%.4f", median(mdd_res$naive_SE_over_donoraware_SE, na.rm=TRUE)))
)
writeLines(report, file.path(out, "GSE102556_donor_audit_report.md"), useBytes=TRUE)

# -------------------------------------------------------------------------
# GSE48350 same-data naive versus donor-aware comparison
# -------------------------------------------------------------------------
ad_scores_all <- read.csv(file.path(root, "05_results", "AD_projection", "AD_projection_module_scores.csv"), check.names=FALSE)
ad_scores <- ad_scores_all[ad_scores_all$dataset == "AD_GSE48350", ]
ad_map <- read.csv(file.path(root, "08_donor_aware_reanalysis", "02_GSE48350_donor_region_map.csv"), check.names=FALSE)
ad <- merge(ad_scores, ad_map[,c("GSM","donor_id","region","diagnosis","age","sex")], by.x="sample_id", by.y="GSM", all.x=TRUE, sort=FALSE,
            suffixes=c("_score",""))
if (nrow(ad) != 253L || anyNA(ad$donor_id)) stop("GSE48350 requires exactly 253 mapped samples")
ad$diagnosis <- factor(ad$diagnosis, levels=c("Control","AD")); ad$region <- factor(ad$region)
ad$sex <- factor(tolower(ad$sex)); ad$donor_id <- factor(ad$donor_id); ad$age_z <- z(as.numeric(ad$age))
analytic_ad <- droplevels(ad[complete.cases(ad[,c("diagnosis","region","age_z","sex","donor_id")]),])
if (nrow(analytic_ad) != 253L) stop("GSE48350 complete-case model did not retain all 253 samples")

fit_ad_module <- function(module) {
  d <- analytic_ad; d$score_z <- z(d[[module]])
  da <- lmer(score_z ~ diagnosis + region + age_z + sex + (1|donor_id), data=d, REML=TRUE, control=lmerControl(optimizer="bobyqa"))
  nv <- lm(score_z ~ diagnosis + region + age_z + sex, data=d)
  eda <- extract_fixed(da,"diagnosisAD"); env <- extract_fixed(nv,"diagnosisAD")
  cr <- tryCatch({
    ct <- clubSandwich::coef_test(nv, vcov="CR2", cluster=d$donor_id, test="Satterthwaite")
    rr <- ct[rownames(ct)=="diagnosisAD",]
    c(beta=rr$beta, SE=rr$SE, df=rr$df_Satt, P=rr$p_Satt)
  }, error=function(e)c(beta=NA,SE=NA,df=NA,P=NA))
  vc <- as.data.frame(VarCorr(da))
  data.frame(module=module,
             naive_beta=env["beta"],naive_SE=env["SE"],naive_df=env["df"],naive_CI_low=env["CI_low"],naive_CI_high=env["CI_high"],naive_P=env["P"],
             donoraware_beta=eda["beta"],donoraware_SE=eda["SE"],donoraware_df=eda["df"],donoraware_CI_low=eda["CI_low"],donoraware_CI_high=eda["CI_high"],donoraware_P=eda["P"],
             cluster_robust_beta=cr["beta"],cluster_robust_SE=cr["SE"],cluster_robust_df=cr["df"],cluster_robust_P=cr["P"],
             sample_n=nrow(d),donor_n=nlevels(d$donor_id),
             random_intercept_variance=vc$vcov[vc$grp=="donor_id"][1],residual_variance=vc$vcov[vc$grp=="Residual"][1],
             singular_fit=isSingular(da,tol=1e-5),
             naive_formula="score_z ~ diagnosis + region + age_z + sex",
             donoraware_formula="score_z ~ diagnosis + region + age_z + sex + (1 | donor_id)",
             cluster_robust_formula="same fixed-effects OLS with donor-clustered CR2 SE",stringsAsFactors=FALSE)
}
ad_res <- do.call(rbind,lapply(module_names,fit_ad_module))
ad_res$naive_BH_FDR <- p.adjust(ad_res$naive_P,"BH")
ad_res$donoraware_BH_FDR <- p.adjust(ad_res$donoraware_P,"BH")
ad_res$cluster_robust_BH_FDR <- p.adjust(ad_res$cluster_robust_P,"BH")
ad_res$beta_difference_donoraware_minus_naive <- ad_res$donoraware_beta-ad_res$naive_beta
ad_res$SE_ratio_donoraware_over_naive <- ad_res$donoraware_SE/ad_res$naive_SE
ad_res$significance_transition <- paste(fdr_status(ad_res$naive_BH_FDR),"to",fdr_status(ad_res$donoraware_BH_FDR))
write_csv(ad_res,"GSE48350_naive_vs_donoraware_full.csv")

# Quantitative-grid figures: direct same-data comparison and uncertainty change.
p_ad_beta <- ggplot(ad_res,aes(naive_beta,donoraware_beta,color=significance_transition))+
  geom_hline(yintercept=0,color="#BDBDBD",linewidth=.3)+geom_vline(xintercept=0,color="#BDBDBD",linewidth=.3)+
  geom_abline(slope=1,intercept=0,linetype=2,color="#777777",linewidth=.35)+geom_point(size=1.7)+
  scale_color_manual(values=c("FDR<0.05 to not_FDR<0.05"="#D24B40","FDR<0.05 to FDR<0.05"="#3182BD","not_FDR<0.05 to FDR<0.05"="#33B5A5","not_FDR<0.05 to not_FDR<0.05"="#767676"),drop=FALSE)+
  coord_equal()+labs(x="Naive diagnosis coefficient (SD)",y="Donor-aware diagnosis coefficient (SD)",color="BH-FDR transition",title="Same-sample effect estimates")
save_plot(p_ad_beta,"Figure_GSE48350_beta_comparison",89,80)

p_ad_se <- ggplot(ad_res,aes(reorder(module,SE_ratio_donoraware_over_naive),SE_ratio_donoraware_over_naive,fill=significance_transition))+
  geom_hline(yintercept=1,linetype=2,color="#777777",linewidth=.35)+geom_col(width=.72)+coord_flip()+
  scale_fill_manual(values=c("FDR<0.05 to not_FDR<0.05"="#D24B40","FDR<0.05 to FDR<0.05"="#3182BD","not_FDR<0.05 to FDR<0.05"="#33B5A5","not_FDR<0.05 to not_FDR<0.05"="#BDBDBD"),drop=FALSE)+
  labs(x=NULL,y="Donor-aware SE / naive SE",fill="BH-FDR transition",title="Uncertainty after donor adjustment")
save_plot(p_ad_se,"Figure_GSE48350_SE_ratio",120,115)

transition_counts <- as.data.frame(table(ad_res$significance_transition)); names(transition_counts)<-c("transition","modules")
p_ad_tr <- ggplot(transition_counts,aes(reorder(transition,modules),modules,fill=transition))+geom_col(width=.65)+coord_flip()+
  scale_fill_manual(values=c("FDR<0.05 to not_FDR<0.05"="#D24B40","FDR<0.05 to FDR<0.05"="#3182BD","not_FDR<0.05 to FDR<0.05"="#33B5A5","not_FDR<0.05 to not_FDR<0.05"="#BDBDBD"),guide="none")+
  labs(x=NULL,y="Modules (n)",title="BH-FDR transitions across 24 modules")
save_plot(p_ad_tr,"Figure_GSE48350_significance_transition",89,70)

# GSE102556 figures
long_mdd <- rbind(
  data.frame(module=mdd_res$module,model="Naive",beta=mdd_res$naive_beta,low=mdd_res$naive_CI_low,high=mdd_res$naive_CI_high,FDR=mdd_res$naive_BH_FDR),
  data.frame(module=mdd_res$module,model="Donor-aware",beta=mdd_res$donoraware_beta,low=mdd_res$donoraware_CI_low,high=mdd_res$donoraware_CI_high,FDR=mdd_res$donoraware_BH_FDR))
ord <- mdd_res$module[order(mdd_res$donoraware_beta)]
long_mdd$module <- factor(long_mdd$module,levels=ord)
p_mdd <- ggplot(long_mdd,aes(beta,module,color=model))+
  geom_vline(xintercept=0,color="#BDBDBD",linewidth=.3)+geom_errorbarh(aes(xmin=low,xmax=high),height=.16,position=position_dodge(width=.5),linewidth=.35)+
  geom_point(position=position_dodge(width=.5),size=1.3)+scale_color_manual(values=c("Naive"="#E28E2C","Donor-aware"="#3182BD"))+
  labs(x="Diagnosis coefficient (SD), 95% CI",y=NULL,color=NULL,title="GSE102556: identical samples and covariates")
save_plot(p_mdd,"Figure_GSE102556_model_comparison",150,120)

region_plot <- region_res
region_plot$module <- factor(region_plot$module,levels=ord)
p_reg <- ggplot(region_plot,aes(beta,module,color=BH_FDR_global_144<.05))+
  geom_vline(xintercept=0,color="#BDBDBD",linewidth=.25)+geom_errorbarh(aes(xmin=CI_low,xmax=CI_high),height=0,linewidth=.25)+geom_point(size=.75)+
  facet_wrap(~brain_region,ncol=3)+scale_color_manual(values=c("FALSE"="#767676","TRUE"="#D24B40"),labels=c("No","Yes"))+
  labs(x="Region-specific diagnosis coefficient (SD), 95% CI",y=NULL,color="Global 144-test BH-FDR<0.05",title="Descriptive MDD brain-region effects")+
  theme(axis.text.y=element_text(size=4.8))
save_plot(p_reg,"Figure_GSE102556_region_effects",183,150)

cat("GSE102556 samples/donors:",nrow(cw),length(unique(cw$donor_id)),"\n")
cat("GSE102556 beta correlation:",cor(mdd_res$naive_beta,mdd_res$donoraware_beta),"\n")
cat("GSE48350 beta correlation:",cor(ad_res$naive_beta,ad_res$donoraware_beta),"\n")
cat("GSE48350 FDR counts naive/donor-aware:",sum(ad_res$naive_BH_FDR<.05),sum(ad_res$donoraware_BH_FDR<.05),"\n")
