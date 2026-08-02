########################################
# Script: 04_preprocessing_audit.R
# Purpose: Audit dataset-specific transformations, missingness, redundant normalization, and ssGSEA sensitivity.
# Input: Staged expression matrices, module definitions, and phenotype tables.
# Output: Preprocessing registry, sensitivity tables, and figures.
# Software: R
# Version: 4.5.2
# Random seed: Not applicable (deterministic)
# Author: Study authors
########################################
options(stringsAsFactors = FALSE, width = 180, contrasts = c("contr.treatment", "contr.poly"))
suppressPackageStartupMessages({
  library(edgeR)
  library(GSVA)
  library(limma)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out <- file.path(root, "13_methodological_reanalysis_v9")
mat_out <- file.path(out, "preprocessed_main_old_mapping")
score_out <- file.path(out, "module_scores_main_old_mapping")
dir.create(mat_out, recursive = TRUE, showWarnings = FALSE)
dir.create(score_out, recursive = TRUE, showWarnings = FALSE)
set.seed(20260802)

write_csv <- function(x, name) utils::write.csv(x, file.path(out, name), row.names = FALSE, na = "")
z <- function(x) as.numeric(scale(x))
num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.+-]", "", as.character(x))))

modules <- readRDS(file.path(root, "04_intermediate", "gene_sets", "liver_modules_combined.rds"))
modules <- lapply(modules, function(x) unique(toupper(trimws(as.character(x)))))

clean_matrix <- function(x) {
  storage.mode(x) <- "double"
  n0s <- ncol(x); n0g <- nrow(x); miss0 <- mean(!is.finite(x))
  sample_missing <- colMeans(!is.finite(x))
  keep_s <- is.finite(sample_missing) & sample_missing <= .20
  x <- x[, keep_s, drop=FALSE]
  keep_g <- rowSums(is.finite(x)) > 0
  x <- x[keep_g,,drop=FALSE]
  imputed <- 0L
  if (anyNA(x)) {
    for (i in which(rowSums(is.na(x)) > 0L)) {
      idx <- is.na(x[i,]); med <- median(x[i,!idx],na.rm=TRUE); x[i,idx] <- med; imputed <- imputed + sum(idx)
    }
  }
  vars <- apply(x,1,var)
  keep_v <- is.finite(vars) & vars > 0
  x <- x[keep_v,,drop=FALSE]
  list(matrix=x, audit=data.frame(initial_samples=n0s,final_samples=ncol(x),removed_samples=n0s-ncol(x),
                                  initial_genes=n0g,final_genes=nrow(x),removed_all_missing_or_zero_variance=n0g-nrow(x),
                                  initial_missing_fraction=miss0,imputed_cells=imputed,imputed_fraction=imputed/(n0s*n0g)))
}

ssgsea_all <- function(expr) {
  shared_sets <- lapply(modules, intersect, y=rownames(expr))
  shared_sets <- shared_sets[lengths(shared_sets)>=3L]
  param <- GSVA::ssgseaParam(as.matrix(expr), shared_sets, minSize=3, normalize=TRUE)
  as.matrix(GSVA::gsva(param, verbose=FALSE, BPPARAM=BiocParallel::SerialParam()))
}

datasets <- data.frame(
  dataset_id=c("liver_GSE135251","liver_GSE167523","liver_GSE126848","MDD_GSE98793","MDD_GSE102556","AD_GSE48350","AD_GSE5281","AD_GSE33000"),
  GEO_accession=c("GSE135251","GSE167523","GSE126848","GSE98793","GSE102556","GSE48350","GSE5281","GSE33000"),
  platform=c("GPL18573","GPL21290","GPL18573","GPL570","GPL11154","GPL570","GPL570","GPL4372"),
  tissue=c("liver","liver","liver","blood","brain (6 regions)","brain (4 regions)","brain (5 regions)","brain (DLPFC)"),
  disease=c("MASLD","MASLD","MASLD","MDD","MDD","AD","AD","AD"),
  technology=c("RNA-seq","RNA-seq","RNA-seq","microarray","RNA-seq","microarray","microarray","custom two-color microarray"),
  input_path=c(
    "04_intermediate/expression_matrices/liver_GSE135251/expression_matrix_raw.rds",
    "04_intermediate/expression_matrices/liver_GSE167523/expression_matrix_raw.rds",
    "04_intermediate/expression_matrices/liver_GSE126848/expression_matrix_raw.rds",
    "04_intermediate/expression_matrices/MDD_GSE98793/expression_matrix_gene_level.rds",
    "04_intermediate/expression_matrices/MDD_GSE102556/expression_matrix_raw.rds",
    "04_intermediate/expression_matrices/AD_GSE48350/expression_matrix_gene_level.rds",
    "04_intermediate/expression_matrices/AD_GSE5281/expression_matrix_gene_level.rds",
    "04_intermediate/expression_matrices/AD_GSE33000/expression_matrix_gene_level.rds"),
  source_type=c("supplementary raw gene counts","supplementary raw gene counts","supplementary raw gene counts","GEO processed matrix","supplementary FPKM matrix","GEO processed matrix","GEO processed matrix","GEO processed Rosetta expression matrix"),
  geo_processing=c(
    "STAR GRCh38; HTSeq-count; supplementary gene-level counts",
    "TopHat/Bowtie2 hg19; Cufflinks FPKM also reported; project input is supplementary raw count matrix",
    "STAR GRCh38; protein-coding filtering; supplementary count table",
    "RMA with Bioconductor affy default parameters",
    "TopHat hg19/GENCODE; HTSeq counts; provided supplementary FPKM matrix",
    "GC-RMA or PLIER; values <0.01 set to 0.01; per-chip 50th percentile and per-gene median normalization",
    "MAS5.0 with target scaling to 150",
    "Rosetta Resolver error modeling and error-weighted squeeze summarization"),
  old_processing=c(
    "log2(count+1) because max>1000","log2(count+1) because max>1000","log2(count+1) because max>1000",
    "additional quantile normalization despite reported RMA","log2(FPKM+1) before generic RNA-seq branch",
    "additional quantile normalization despite reported GC-RMA/PLIER normalization",
    "additional quantile normalization despite MAS5 target scaling","additional quantile normalization of Rosetta-processed values"),
  audited_main=c(
    "TMM library-size normalization followed by edgeR logCPM (prior.count=1)",
    "TMM library-size normalization followed by edgeR logCPM (prior.count=1)",
    "TMM library-size normalization followed by edgeR logCPM (prior.count=1)",
    "retain GEO RMA processed values; no second quantile normalization",
    "log2(FPKM+1); no quantile normalization",
    "retain GEO GC-RMA/PLIER processed values; no second quantile normalization",
    "retain GEO MAS5 target-scaled processed values; no second quantile normalization",
    "retain GEO Rosetta processed values; no second quantile normalization"),
  stringsAsFactors=FALSE
)

matrix_inventory <- read.csv(file.path(out,"input_expression_matrix_inventory.csv"),check.names=FALSE)
audit_rows <- list(); registry_rows <- list(); score_paths <- character()
for(i in seq_len(nrow(datasets))){
  d <- datasets[i,]
  raw <- readRDS(file.path(root,d$input_path))
  clean <- clean_matrix(raw); x <- clean$matrix
  if(grepl("liver_GSE",d$dataset_id)){
    y <- edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(counts=x),method="TMM"),log=TRUE,prior.count=1)
  } else if(d$dataset_id=="MDD_GSE102556"){
    y <- log2(x+1)
  } else {
    y <- x
  }
  saveRDS(y,file.path(mat_out,paste0(d$dataset_id,".rds")))
  sc <- ssgsea_all(y)
  saveRDS(sc,file.path(score_out,paste0(d$dataset_id,".rds")))
  score_paths[d$dataset_id] <- file.path(score_out,paste0(d$dataset_id,".rds"))
  vals<-as.numeric(raw); finite<-vals[is.finite(vals)]; qs<-quantile(finite,c(0,.25,.5,.75,1),names=FALSE)
  audit_rows[[i]]<-cbind(data.frame(dataset_id=d$dataset_id),clean$audit)
  registry_rows[[i]]<-data.frame(
    GEO_accession=d$GEO_accession,GPL=d$platform,tissue=d$tissue,disease=d$disease,
    download_object=ifelse(grepl("GEO processed",d$source_type),"GEO series/project-cached processed matrix","downloaded supplementary expression table"),
    expression_matrix_source=d$input_path,raw_or_processed=d$source_type,technology=d$technology,
    original_paper_preprocessing=d$geo_processing,GEO_metadata_preprocessing=d$geo_processing,
    min=qs[1],Q1=qs[2],median=qs[3],Q3=qs[4],max=qs[5],has_negative=any(finite<0),
    suspected_log2_scale=ifelse(d$dataset_id %in% c("MDD_GSE98793","AD_GSE33000"),"yes",ifelse(grepl("liver|102556",d$dataset_id),"no (count/FPKM input)","no; positive processed scale")),
    suspected_already_quantile_normalized=ifelse(d$dataset_id %in% c("MDD_GSE98793","AD_GSE48350"),"yes/report indicates normalized",ifelse(d$dataset_id=="AD_GSE5281","MAS5 target-scaled, not assumed quantile-normalized",ifelse(d$dataset_id=="AD_GSE33000","processed by Rosetta; no evidence supporting another quantile normalization","not applicable"))),
    old_code_transformation=d$old_processing,audited_main_processing=d$audited_main,
    additional_quantile_normalization_in_main=FALSE,final_sample_n=ncol(y),final_gene_n=nrow(y),
    removed_samples_and_reason=ifelse(clean$audit$removed_samples==0,"0","samples with >20% missing values"),
    missing_fraction=clean$audit$initial_missing_fraction,missing_value_handling=paste0("feature-median imputation: ",clean$audit$imputed_cells," cells (",signif(clean$audit$imputed_fraction,3),")"),
    mapping_strategy=ifelse(grepl("microarray",tolower(d$technology)),"legacy first-symbol mapping for Phase 3; strict mapping evaluated separately in Phase 4","gene symbols supplied by expression table; duplicate symbols collapsed by prespecified feature rule where needed"),
    notes="No raw file was modified; all audited matrices were written to the versioned v9 directory.",stringsAsFactors=FALSE)
  rm(raw,x,y,sc,vals,finite);invisible(gc())
}
write_csv(do.call(rbind,audit_rows),"missing_value_audit.csv")

# Add GSE144136 marker-resource row (not a full expression-matrix projection).
registry_rows[[length(registry_rows)+1L]]<-data.frame(
  GEO_accession="GSE144136",GPL="GPL20301",tissue="brain DLPFC",disease="MDD marker resource",
  download_object="published supplementary between-cluster marker table and project metadata",expression_matrix_source="01_raw_data/single_cell_snRNA/GSE144136/expression_matrix.rds (0-feature placeholder); marker analysis used external supplementary table",
  raw_or_processed="published processed marker result",technology="single-nucleus RNA-seq",
  original_paper_preprocessing="10x single-nucleus RNA-seq; published cluster marker analysis used as an annotation resource",
  GEO_metadata_preprocessing="full count matrix not represented by the 0-feature project placeholder",
  min=NA,Q1=NA,median=NA,Q3=NA,max=NA,has_negative=NA,suspected_log2_scale="not evaluable",suspected_already_quantile_normalized="not applicable",
  old_code_transformation="no module-score expression transformation; marker table filters only",audited_main_processing="retain role as published marker resource; do not treat 34 GEO entries as a full expression matrix or nuclei",
  additional_quantile_normalization_in_main=FALSE,final_sample_n=34,final_gene_n=NA,removed_samples_and_reason="not applicable",missing_fraction=NA,
  missing_value_handling="not applicable",mapping_strategy="reported gene symbols in supplementary marker table",notes="Not used in AD/MDD module-effect modeling.",stringsAsFactors=FALSE)
registry <- do.call(rbind,registry_rows)
write_csv(registry,"dataset_preprocessing_registry.csv")

# -------------------------------------------------------------------------
# Additional quantile-normalization sensitivity for all four microarrays.
# -------------------------------------------------------------------------
micro_ids <- c("MDD_GSE98793","AD_GSE48350","AD_GSE5281","AD_GSE33000")
pheno_path <- function(id) file.path(root,"04_intermediate","phenotype_tables",id,"pheno_harmonized.rds")
ad_map <- read.csv(file.path(root,"08_donor_aware_reanalysis","02_GSE48350_donor_region_map.csv"),check.names=FALSE)

extract <- function(fit,term){s<-coef(summary(fit));pcol<-grep("Pr\\(",colnames(s),value=TRUE)[1];c(beta=s[term,"Estimate"],SE=s[term,"Std. Error"],P=s[term,pcol])}
fit_effects <- function(id,score_mat){
  ph<-readRDS(pheno_path(id)); common<-intersect(colnames(score_mat),ph$sample_id); score_mat<-score_mat[,common,drop=FALSE];ph<-ph[match(common,ph$sample_id),,drop=FALSE]
  outl<-lapply(rownames(score_mat),function(m){
    d<-data.frame(score_z=z(score_mat[m,]),stringsAsFactors=FALSE)
    if(id=="AD_GSE48350"){
      mp<-ad_map[match(common,ad_map$GSM),];d$diagnosis<-factor(mp$diagnosis,levels=c("Control","AD"));d$region<-factor(mp$region);d$age_z<-z(mp$age);d$sex<-factor(tolower(mp$sex));d$donor_id<-factor(mp$donor_id)
      fit<-lmer(score_z~diagnosis+region+age_z+sex+(1|donor_id),data=d,REML=TRUE,control=lmerControl(optimizer="bobyqa"));e<-extract(fit,"diagnosisAD");formula<-"score_z ~ diagnosis + region + age_z + sex + (1 | donor_id)"
    } else if(id=="AD_GSE33000"){
      d$diagnosis<-factor(ph$diagnosis,levels=c("Control","AD"));d$age_z<-z(num(ph[["age:ch2"]]));d$sex<-factor(tolower(ph[["gender:ch2"]]));d<-d[!is.na(d$diagnosis),];fit<-lm(score_z~diagnosis+age_z+sex,data=d);e<-extract(fit,"diagnosisAD");formula<-"score_z ~ diagnosis + age_z + sex"
    } else if(id=="AD_GSE5281"){
      d$diagnosis<-factor(ph$diagnosis,levels=c("Control","AD"));d$region<-factor(ph$brain_region);fit<-lm(score_z~diagnosis+region,data=d);e<-extract(fit,"diagnosisAD");formula<-"descriptive score_z ~ diagnosis + region"
    } else {
      d$diagnosis<-factor(ph$diagnosis,levels=c("Control","MDD"));fit<-lm(score_z~diagnosis,data=d);e<-extract(fit,"diagnosisMDD");formula<-"score_z ~ diagnosis"
    }
    data.frame(dataset_id=id,module=m,beta=e["beta"],SE=e["SE"],P=e["P"],sample_n=nrow(d),formula=formula)
  });do.call(rbind,outl)
}

dn_rows<-list(); eff_rows<-list()
for(id in micro_ids){
  x<-readRDS(file.path(mat_out,paste0(id,".rds")))
  xq<-limma::normalizeBetweenArrays(x,method="quantile")
  s0<-readRDS(file.path(score_out,paste0(id,".rds")));sq<-ssgsea_all(xq)
  e0<-fit_effects(id,s0);eq<-fit_effects(id,sq);e0$BH_FDR<-p.adjust(e0$P,"BH");eq$BH_FDR<-p.adjust(eq$P,"BH")
  mm<-merge(e0,eq,by=c("dataset_id","module"),suffixes=c("_processed_main","_additional_quantile"))
  mm$score_correlation<-vapply(mm$module,function(m)cor(s0[m,],sq[m,],use="pairwise.complete.obs"),numeric(1))
  mm$FDR_status_change<-(mm$BH_FDR_processed_main<.05)!=(mm$BH_FDR_additional_quantile<.05)
  dn_rows[[id]]<-mm
  eff_rows[[id]]<-data.frame(dataset_id=id,mean_module_score_correlation=mean(mm$score_correlation),min_module_score_correlation=min(mm$score_correlation),
                             beta_correlation=cor(mm$beta_processed_main,mm$beta_additional_quantile),
                             median_SE_ratio_additional_over_main=median(mm$SE_additional_quantile/mm$SE_processed_main),
                             significant_modules_main=sum(mm$BH_FDR_processed_main<.05),significant_modules_additional_quantile=sum(mm$BH_FDR_additional_quantile<.05),
                             FDR_status_changes=sum(mm$FDR_status_change),conclusion=ifelse(cor(mm$beta_processed_main,mm$beta_additional_quantile)>.95,"effect estimates broadly stable","material preprocessing sensitivity"))
  rm(x,xq,s0,sq,e0,eq,mm);invisible(gc())
}
double_audit<-do.call(rbind,dn_rows);summary_audit<-do.call(rbind,eff_rows)
write_csv(double_audit,"double_normalization_audit.csv")
write_csv(summary_audit,"preprocessing_sensitivity_summary.csv")

pd<-double_audit
p<-ggplot(pd,aes(beta_processed_main,beta_additional_quantile,color=dataset_id))+
  geom_hline(yintercept=0,color="#BDBDBD",linewidth=.3)+geom_vline(xintercept=0,color="#BDBDBD",linewidth=.3)+geom_abline(slope=1,intercept=0,linetype=2,color="#777777",linewidth=.35)+geom_point(size=1.2,alpha=.8)+
  facet_wrap(~dataset_id,scales="free")+scale_color_manual(values=c("AD_GSE33000"="#3182BD","AD_GSE48350"="#33B5A5","AD_GSE5281"="#E28E2C","MDD_GSE98793"="#D24B40"),guide="none")+
  theme_classic(base_size=7,base_family="Arial")+labs(x="Audited processed-value effect (SD)",y="Additional-quantile effect (SD)",title="Preprocessing sensitivity across microarray datasets")
ggsave(file.path(out,"Figure_preprocessing_sensitivity.pdf"),p,width=183/25.4,height=120/25.4,device=grDevices::cairo_pdf,family="Arial")
ggsave(file.path(out,"Figure_preprocessing_sensitivity.png"),p,width=183/25.4,height=120/25.4,dpi=600,device=ragg::agg_png)

cat("Wrote registry rows:",nrow(registry),"\n")
print(summary_audit)
