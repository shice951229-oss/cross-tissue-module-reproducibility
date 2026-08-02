########################################
# Script: 09_matched_random_modules_MDD.R
# Purpose: Run separate MDD blood and donor-aware brain matched-random-module controls.
# Input: Fixed modules, liver reference, MDD matrices, phenotypes, and donor crosswalk.
# Output: MDD empirical module and global negative-control results.
# Software: R
# Version: 4.5.2
# Random seed: 20260802; B=2000 per module and tissue
# Author: Study authors
########################################
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 180, contrasts = c("contr.treatment", "contr.poly"))
suppressPackageStartupMessages({
  library(GSVA)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out <- file.path(root, "13_methodological_reanalysis_v9")
B <- as.integer(Sys.getenv("MDD_NEG_CONTROL_B", "2000"))
if (!is.finite(B) || B < 1L) stop("Invalid MDD_NEG_CONTROL_B")
seed <- 20260802L
set.seed(seed)
prefix <- if (B >= 2000L) "" else paste0("benchmark_B", B, "_")
z <- function(x) as.numeric(scale(x))

modules <- readRDS(file.path(root, "04_intermediate", "gene_sets", "liver_modules_combined.rds"))
modules <- lapply(modules, function(x) unique(toupper(trimws(as.character(x)))))
ref <- readRDS(file.path(out, "preprocessed_main_old_mapping", "liver_GSE135251.rds"))
targets <- list(
  GSE98793_blood = readRDS(file.path(out, "preprocessed_strict_mapping", "MDD_GSE98793.rds")),
  GSE102556_brain = readRDS(file.path(out, "preprocessed_main_old_mapping", "MDD_GSE102556.rds"))
)

make_sets <- function(target) {
  universe <- intersect(rownames(ref), rownames(target))
  r <- ref[universe, , drop = FALSE]
  mean_ref <- rowMeans(r)
  iqr_ref <- apply(r, 1, IQR)
  dec <- function(v) pmin(10L, pmax(1L, ceiling(rank(v, ties.method = "average") / length(v) * 10)))
  bin <- paste0("M", dec(mean_ref), "_V", dec(iqr_ref)); names(bin) <- universe
  pool <- split(universe, bin)
  sets <- list(); design <- list()
  for (m in names(modules)) {
    g <- intersect(modules[[m]], universe)
    prof <- table(bin[g])
    sets[[paste0(m, "__observed")]] <- g
    for (b in seq_len(B)) {
      sets[[paste0(m, "__null_", b)]] <- unlist(lapply(names(prof), function(k) {
        sample(pool[[k]], as.integer(prof[[k]]), replace = FALSE)
      }), use.names = FALSE)
    }
    design[[m]] <- data.frame(module = m, original_size = length(modules[[m]]),
                              detectable_size = length(g), coverage_fraction = length(g) / length(modules[[m]]),
                              B = B, seed = seed)
  }
  list(universe = universe, sets = sets, design = do.call(rbind, design))
}

score_sets <- function(expr, sets) {
  p <- GSVA::ssgseaParam(as.matrix(expr), sets, minSize = 3, normalize = TRUE)
  as.matrix(GSVA::gsva(p, verbose = FALSE, BPPARAM = BiocParallel::SerialParam()))
}

fit_blood <- function(scores) {
  ph <- readRDS(file.path(root, "04_intermediate", "phenotype_tables", "MDD_GSE98793", "pheno_harmonized.rds"))
  ph <- ph[match(colnames(scores), ph$sample_id), ]
  i1 <- ph$diagnosis == "MDD"; i0 <- ph$diagnosis == "Control"
  n1 <- sum(i1); n0 <- sum(i0)
  m1 <- rowMeans(scores[, i1, drop = FALSE]); m0 <- rowMeans(scores[, i0, drop = FALSE])
  v1 <- apply(scores[, i1, drop = FALSE], 1, var); v0 <- apply(scores[, i0, drop = FALSE], 1, var)
  sp <- sqrt(((n1 - 1) * v1 + (n0 - 1) * v0) / (n1 + n0 - 2))
  effect <- (m1 - m0) / sp
  se <- sqrt(1 / n1 + 1 / n0 + effect^2 / (2 * (n1 + n0)))
  welch_se <- sqrt(v1 / n1 + v0 / n0)
  df <- (v1 / n1 + v0 / n0)^2 / ((v1 / n1)^2 / (n1 - 1) + (v0 / n0)^2 / (n0 - 1))
  P <- 2 * pt(-abs((m1 - m0) / welch_se), df = df)
  data.frame(effect = effect, SE = se, P = P, singular = NA)
}

brain_data <- function(scores) {
  cw <- read.csv(file.path(out, "GSE102556_sample_donor_crosswalk.csv"), check.names = FALSE)
  cw <- cw[match(colnames(scores), cw$sample_id), ]
  data.frame(
    diagnosis = factor(cw$diagnosis, levels = c("Control", "MDD")),
    brain_region = factor(cw$brain_region),
    age_z = z(cw$age), sex = factor(tolower(cw$sex)), RIN_z = z(cw$RIN), PMI_z = z(cw$PMI),
    donor_id = factor(cw$donor_id)
  )
}

fit_brain <- function(scores) {
  d <- brain_data(scores)
  fit0 <- lmerTest::lmer(z(scores[1, ]) ~ diagnosis + brain_region + age_z + sex + RIN_z + PMI_z + (1 | donor_id),
                        data = d, REML = TRUE, control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
  rows <- lapply(seq_len(nrow(scores)), function(i) {
    fit <- if (i == 1L) fit0 else suppressMessages(lme4::refit(fit0, newresp = z(scores[i, ])))
    ss <- coef(summary(fit)); term <- "diagnosisMDD"
    b <- ss[term, "Estimate"]; se <- ss[term, "Std. Error"]
    p <- if ("Pr(>|t|)" %in% colnames(ss)) ss[term, "Pr(>|t|)"] else 2 * pnorm(-abs(b / se))
    c(effect = b, SE = se, P = p, singular = isSingular(fit, tol = 1e-5))
  })
  x <- as.data.frame(do.call(rbind, rows)); x[] <- lapply(x, as.numeric); x
}

all_results <- list(); all_null <- list(); global_rows <- list(); design_rows <- list()
t0 <- Sys.time()
for (dataset in names(targets)) {
  set.seed(seed)
  target <- targets[[dataset]]
  spec <- make_sets(target)
  design <- spec$design; design$dataset <- dataset; design$universe_size <- length(spec$universe)
  design_rows[[dataset]] <- design
  message("Scoring ", length(spec$sets), " sets for ", dataset, "; universe=", length(spec$universe))
  scores <- score_sets(target[spec$universe, , drop = FALSE], spec$sets)
  fits <- if (dataset == "GSE98793_blood") fit_blood(scores) else fit_brain(scores)

  obs_rows <- list(); null_rows <- list()
  for (m in names(modules)) {
    ids <- c(match(paste0(m, "__observed"), rownames(scores)),
             match(paste0(m, "__null_", seq_len(B)), rownames(scores)))
    met <- fits[ids, , drop = FALSE]
    obs <- met[1, ]; nul <- met[-1, ]
    emp <- (1 + sum(abs(nul$effect) >= abs(obs$effect), na.rm = TRUE)) / (B + 1)
    dd <- design[design$module == m, ]
    obs_rows[[m]] <- data.frame(
      dataset = dataset, tissue = ifelse(grepl("blood", dataset), "blood", "brain"), module = m,
      original_size = dd$original_size, detectable_size = dd$detectable_size, coverage_fraction = dd$coverage_fraction,
      effect = obs$effect, SE = obs$SE, CI_low = obs$effect - 1.96 * obs$SE,
      CI_high = obs$effect + 1.96 * obs$SE, P = obs$P,
      empirical_P_abs_effect = emp, singular_fit = as.logical(obs$singular), B = B, seed = seed,
      formula = ifelse(dataset == "GSE98793_blood", "Welch two-sample t-test; Cohen d = MDD minus Control",
                       "score_z ~ diagnosis + brain_region + age_z + sex + RIN_z + PMI_z + (1 | donor_id)")
    )
    null_rows[[m]] <- data.frame(dataset = dataset, module = m, replicate = seq_len(B),
                                 effect = nul$effect, P = nul$P)
  }
  rr <- do.call(rbind, obs_rows); nn <- do.call(rbind, null_rows)
  rr$empirical_BH_FDR_abs_effect <- p.adjust(rr$empirical_P_abs_effect, "BH")
  all_results[[dataset]] <- rr; all_null[[dataset]] <- nn

  glob <- do.call(rbind, lapply(seq_len(B), function(b) {
    x <- nn[nn$replicate == b, ]
    data.frame(replicate = b, mean_abs_effect = mean(abs(x$effect)), nominal_support_count = sum(x$P < .05, na.rm = TRUE))
  }))
  obs_mean <- mean(abs(rr$effect)); obs_count <- sum(rr$P < .05, na.rm = TRUE)
  global_rows[[dataset]] <- data.frame(
    analysis = dataset, metric = c("mean_abs_effect", "nominal_support_count"),
    observed = c(obs_mean, obs_count), null_mean = c(mean(glob$mean_abs_effect), mean(glob$nominal_support_count)),
    null_q025 = c(quantile(glob$mean_abs_effect, .025), quantile(glob$nominal_support_count, .025)),
    null_q975 = c(quantile(glob$mean_abs_effect, .975), quantile(glob$nominal_support_count, .975)),
    tail = "higher",
    empirical_P = c((1 + sum(glob$mean_abs_effect >= obs_mean)) / (B + 1),
                    (1 + sum(glob$nominal_support_count >= obs_count)) / (B + 1))
  )
  message("Completed ", dataset, " after ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
}

res <- do.call(rbind, all_results); null <- do.call(rbind, all_null); glob <- do.call(rbind, global_rows)
glob$BH_FDR <- ave(glob$empirical_P, glob$analysis, FUN = function(p) p.adjust(p, "BH"))
write.csv(res, file.path(out, paste0(prefix, "matched_random_module_results_MDD.csv")), row.names = FALSE, na = "")
write.csv(do.call(rbind, design_rows), file.path(out, paste0(prefix, "matched_random_module_design_MDD_table.csv")), row.names = FALSE, na = "")
write.csv(glob, file.path(out, paste0(prefix, "global_negative_control_tests_MDD.csv")), row.names = FALSE, na = "")
saveRDS(null, file.path(out, paste0(prefix, "matched_random_module_null_MDD.rds")))

writeLines(c(
  "# Matched random-module negative-control design: MDD", "",
  paste0("- Fixed seed: ", seed, "; B=", B, " random sets per observed module and tissue."),
  "- Matching was conducted separately for MDD blood and MDD brain using module size, GSE135251 mean-expression decile, GSE135251 IQR decile, and target detectability.",
  "- GSE98793 used the final strict unambiguous/IQR probe mapping and the prespecified Welch/Cohen-d analysis.",
  "- GSE102556 used the verified 48-donor structure and the donor-aware mixed model.",
  "- Blood and brain were tested and reported separately; no cross-tissue pooled MDD estimate was calculated.",
  "- Empirical P values use (1 + at-least-as-extreme null count)/(B + 1); BH correction spans 24 modules within each tissue."
), file.path(out, paste0(prefix, "matched_random_module_design_MDD.md")))

if (B >= 2000L) {
  pd <- merge(null, res[, c("dataset", "module", "effect")], by = c("dataset", "module"), suffixes = c("_null", "_observed"))
  pd$module <- factor(pd$module, levels = rev(unique(res$module)))
  p <- ggplot(pd, aes(abs(effect_null), module)) +
    geom_violin(fill = "#D8D8D8", color = NA, scale = "width") +
    geom_point(data = res, aes(abs(effect), module), inherit.aes = FALSE, color = "#D24B40", size = 1.1) +
    facet_wrap(~dataset, scales = "free_x") + theme_classic(base_size = 6.5, base_family = "Arial") +
    labs(x = "Absolute MDD effect under matched null", y = NULL,
         title = "MDD projections against tissue-specific matched random modules")
  ggsave(file.path(out, "Figure_random_null_MDD_effects.pdf"), p, width = 180 / 25.4, height = 125 / 25.4,
         device = grDevices::cairo_pdf, family = "Arial")
  ggsave(file.path(out, "Figure_random_null_MDD_effects.png"), p, width = 180 / 25.4, height = 125 / 25.4,
         dpi = 600, device = ragg::agg_png)
}

message("Completed MDD matched-random controls with B=", B)
