########################################
# Script: 17_update_figures_v10.R
# Purpose: Apply text-only refinements to Figures 1 and 6 from locked tables.
# Input: Integrated evidence and liver-preservation CSV files.
# Output: Final Figure 1 and Figure 6 PDF/PNG files.
# Software: R
# Version: 4.5.2
# Random seed: Not applicable (deterministic)
# Author: Study authors
########################################
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

root <- normalizePath(".", winslash = "/")
out <- file.path(root, "13_methodological_reanalysis_v9")

save_plot <- function(p, stem, width_mm, height_mm) {
  ggsave(file.path(out, paste0(stem, ".pdf")), p, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(file.path(out, paste0(stem, ".png")), p, width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
}

pretty_module <- function(x) {
  x <- sub("^hepatic_", "", x)
  x <- sub("_module$", "", x)
  x <- gsub("_", " ", x)
  x <- sub("^WGCNA ", "WGCNA ", x)
  x
}

# Figure 1: preserve the v9 layout and replace only the inferential-boundary wording.
nodes <- data.frame(
  x = c(1, 1, 1, 2.6, 4.2, 5.8),
  y = c(3.3, 2.2, 1.1, 2.2, 2.2, 2.2),
  label = c(
    "Liver reference\nGSE135251",
    "External liver\nGSE126848 / GSE167523",
    "HPA tissue\nconsensus",
    "24 fixed modules\nssGSEA within dataset",
    "AD primary + MDD external\ndonor-aware / naive sensitivity",
    "Matched-random modules\nEvidence tiers + limits"
  ),
  group = c("input", "input", "input", "score", "model", "inference")
)
edges <- data.frame(
  x = c(1.45, 1.45, 1.45, 3.05, 4.65),
  y = c(3.3, 2.2, 1.1, 2.2, 2.2),
  xend = c(2.15, 2.15, 2.15, 3.75, 5.35),
  yend = c(2.45, 2.2, 1.95, 2.2, 2.2)
)
p1 <- ggplot() +
  geom_segment(data = edges, aes(x, y, xend = xend, yend = yend),
               arrow = arrow(length = grid::unit(2, "mm")), color = "#767676") +
  geom_label(data = nodes, aes(x, y, label = label, fill = group), size = 3,
             linewidth = .3, label.padding = grid::unit(2.5, "mm"), lineheight = .95) +
  scale_fill_manual(values = c(input = "#EAF2F8", score = "#E8F5E9", model = "#FFF3E0", inference = "#FDEDEC")) +
  annotate("text", x = 3.4, y = 0.60,
           label = "Cross-tissue association does not establish tissue origin,\ninter-organ transfer, formal mediation, or causality.",
           size = 3.15, lineheight = 1.05, color = "#7F2704") +
  coord_cartesian(xlim = c(.35, 6.45), ylim = c(.25, 3.75), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(legend.position = "none", plot.title = element_text(face = "bold", size = 12)) +
  labs(title = "Audit-first cross-tissue module projection and evidence hierarchy")
save_plot(p1, "Figure1_study_design_evidence_hierarchy_v10", 180, 100)

# Figure 6: redraw from locked v9 evidence and preservation tables; no tier or statistic is recalculated.
ev <- read.csv(file.path(out, "integrated_module_evidence_matrix.csv"), check.names = FALSE)
pres <- read.csv(file.path(out, "liver_module_preservation_results.csv"), check.names = FALSE)
tier_levels <- c("Tier D", "Tier C", "Tier B", "Tier A")
ev$final_evidence_tier <- factor(ev$final_evidence_tier, levels = tier_levels)
ev$module_label <- pretty_module(ev$module)
ev$module_label <- factor(ev$module_label, levels = ev$module_label[order(ev$final_evidence_tier)])
p6a <- ggplot(ev, aes(final_evidence_tier, module_label, fill = final_evidence_tier)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("Tier A" = "#1B7837", "Tier B" = "#5AAE61", "Tier C" = "#FDB863", "Tier D" = "#BDBDBD")) +
  theme_minimal(base_size = 6.5, base_family = "Arial") +
  theme(legend.position = "none", panel.grid = element_blank()) +
  labs(x = NULL, y = NULL, title = "A  Final conservative tier")

pres$module_label <- pretty_module(pres$module)
p6b <- ggplot(pres, aes(Zsummary, module_label, color = preservation_class, size = common_gene_n)) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "#767676") +
  geom_vline(xintercept = 10, linetype = "dotted", color = "#4D4D4D") +
  geom_point(alpha = .85) +
  facet_wrap(~test_dataset) +
  scale_color_manual(values = c("no clear evidence" = "#BDBDBD", moderate = "#3182BD", strong = "#D24B40", "low common-gene count" = "#7B3294")) +
  scale_size_continuous(range = c(1, 4)) +
  guides(size = guide_legend(order = 1, nrow = 1), color = guide_legend(order = 2, nrow = 2, byrow = TRUE)) +
  theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 5.5)) +
  labs(x = "Preservation Zsummary", y = NULL, color = NULL, size = "Common genes", title = "B  External liver-network preservation")

p6 <- p6a + p6b + plot_layout(widths = c(.55, 1.45)) +
  plot_annotation(title = "Evidence hierarchy and liver-related attribute checks")
save_plot(p6, "Figure6_final_evidence_and_liver_preservation_v10", 180, 155)

message("V10 Figure 1 and Figure 6 text-only refinements completed from locked result tables.")
