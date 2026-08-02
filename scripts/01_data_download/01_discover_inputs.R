########################################
# Script: 01_discover_inputs.R
# Purpose: Discover and summarize locked expression and phenotype inputs without modifying them.
# Input: Project-staged RDS inputs and phenotype tables.
# Output: Input and phenotype inventories.
# Software: R
# Version: 4.5.2
# Random seed: Not applicable (deterministic)
# Author: Study authors
########################################
options(stringsAsFactors = FALSE, width = 200)

root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "13_methodological_reanalysis_v9")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")

as_numeric_matrix <- function(x) {
  if (inherits(x, "ExpressionSet")) return(Biobase::exprs(x))
  if (is.matrix(x)) return(x)
  if (is.data.frame(x)) return(as.matrix(x))
  if (is.list(x)) {
    candidates <- c("expr", "expression", "expression_matrix", "matrix", "counts", "data")
    for (nm in intersect(candidates, names(x))) {
      y <- x[[nm]]
      if (is.matrix(y) || is.data.frame(y)) return(as.matrix(y))
    }
  }
  NULL
}

matrix_files <- list.files(root, pattern = "expression_matrix.*\\.rds$", recursive = TRUE, full.names = TRUE)
matrix_files <- matrix_files[grepl("^(01_raw_data|04_intermediate)/", sub(paste0("^", gsub("([\\W])", "\\\\\\1", root), "/?"), "", gsub("\\\\", "/", matrix_files)))]

matrix_rows <- vector("list", length(matrix_files))
for (i in seq_along(matrix_files)) {
  f <- matrix_files[[i]]
  rel <- substring(normalizePath(f, winslash = "/"), nchar(root) + 2L)
  obj <- tryCatch(readRDS(f), error = identity)
  if (inherits(obj, "error")) {
    matrix_rows[[i]] <- data.frame(relative_path = rel, read_status = "error", error = conditionMessage(obj))
    next
  }
  mat <- as_numeric_matrix(obj)
  if (is.null(mat)) {
    matrix_rows[[i]] <- data.frame(relative_path = rel, read_status = "not_matrix_like", object_class = paste(class(obj), collapse = ";"), object_names = paste(names(obj), collapse = ";"))
    next
  }
  storage.mode(mat) <- "double"
  vals <- as.numeric(mat)
  finite <- vals[is.finite(vals)]
  qs <- if (length(finite)) unname(stats::quantile(finite, c(0, .25, .5, .75, 1), na.rm = TRUE, names = FALSE)) else rep(NA_real_, 5)
  sample_vals <- if (length(finite) > 1000000L) finite[seq.int(1L, length(finite), length.out = 1000000L)] else finite
  integer_fraction <- if (length(sample_vals)) mean(abs(sample_vals - round(sample_vals)) < 1e-8) else NA_real_
  row_ids <- rownames(mat); col_ids <- colnames(mat)
  matrix_rows[[i]] <- data.frame(
    relative_path = rel,
    read_status = "ok",
    object_class = paste(class(obj), collapse = ";"),
    matrix_class = paste(class(mat), collapse = ";"),
    n_features = nrow(mat), n_samples = ncol(mat),
    min = qs[1], q1 = qs[2], median = qs[3], q3 = qs[4], max = qs[5],
    has_negative = any(finite < 0),
    missing_fraction = mean(!is.finite(vals)),
    integer_like_fraction_sampled = integer_fraction,
    row_id_example = paste(utils::head(row_ids, 3), collapse = ";"),
    col_id_example = paste(utils::head(col_ids, 3), collapse = ";"),
    duplicated_row_ids = if (is.null(row_ids)) NA_integer_ else sum(duplicated(row_ids)),
    duplicated_col_ids = if (is.null(col_ids)) NA_integer_ else sum(duplicated(col_ids)),
    stringsAsFactors = FALSE
  )
  rm(obj, mat, vals, finite, sample_vals); invisible(gc())
}
all_names <- unique(unlist(lapply(matrix_rows, names)))
matrix_rows <- lapply(matrix_rows, function(d) { for (nm in setdiff(all_names, names(d))) d[[nm]] <- NA; d[all_names] })
write_csv(do.call(rbind, matrix_rows), file.path(out_dir, "input_expression_matrix_inventory.csv"))

pheno_files <- list.files(file.path(root, "04_intermediate", "phenotype_tables"), pattern = "pheno_harmonized\\.rds$", recursive = TRUE, full.names = TRUE)
pheno_rows <- list(); dataset_rows <- list()
for (f in pheno_files) {
  rel <- substring(normalizePath(f, winslash = "/"), nchar(root) + 2L)
  dataset_id <- basename(dirname(f))
  x <- readRDS(f)
  if (inherits(x, "AnnotatedDataFrame")) x <- Biobase::pData(x)
  if (!is.data.frame(x)) x <- as.data.frame(x)
  dataset_rows[[length(dataset_rows) + 1L]] <- data.frame(
    dataset_id = dataset_id, relative_path = rel, n_rows = nrow(x), n_columns = ncol(x),
    rowname_example = paste(head(rownames(x), 3), collapse = ";"), stringsAsFactors = FALSE
  )
  for (nm in names(x)) {
    v <- x[[nm]]
    examples <- unique(as.character(v[!is.na(v) & nzchar(as.character(v))]))
    pheno_rows[[length(pheno_rows) + 1L]] <- data.frame(
      dataset_id = dataset_id, column = nm, class = paste(class(v), collapse = ";"),
      missing_fraction = mean(is.na(v) | !nzchar(trimws(as.character(v)))),
      n_unique_nonmissing = length(examples), examples = paste(head(examples, 5), collapse = " | "),
      stringsAsFactors = FALSE
    )
  }
}
write_csv(do.call(rbind, dataset_rows), file.path(out_dir, "phenotype_table_inventory.csv"))
write_csv(do.call(rbind, pheno_rows), file.path(out_dir, "phenotype_column_inventory.csv"))

modules_path <- file.path(root, "04_intermediate", "gene_sets", "liver_modules_combined.rds")
modules <- readRDS(modules_path)
if (!is.list(modules)) stop("Combined module object is not a list")
module_rows <- do.call(rbind, lapply(names(modules), function(nm) {
  genes <- unique(toupper(trimws(as.character(modules[[nm]]))))
  genes <- genes[nzchar(genes)]
  data.frame(module = nm, gene = genes, stringsAsFactors = FALSE)
}))
module_summary <- aggregate(gene ~ module, module_rows, length)
names(module_summary)[2] <- "n_genes"
write_csv(module_rows, file.path(out_dir, "module_gene_lists_24_long.csv"))
write_csv(module_summary, file.path(out_dir, "module_gene_lists_24_summary.csv"))

cat("Expression files:", length(matrix_files), "\n")
cat("Phenotype files:", length(pheno_files), "\n")
cat("Modules:", length(modules), "\n")
