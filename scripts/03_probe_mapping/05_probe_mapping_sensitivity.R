########################################
# Script: 05_probe_mapping_sensitivity.R
# Purpose: Compare legacy and strict unambiguous probe-to-gene mappings.
# Input: Microarray expression matrices, platform annotation, and fixed modules.
# Output: Probe audit, gene coverage, model sensitivity, and comparison figure.
# Software: R
# Version: 4.5.2
# Random seed: Not applicable (deterministic)
# Author: Study authors
########################################
options(stringsAsFactors=FALSE,width=180,contrasts=c("contr.treatment","contr.poly"))
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(org.Hs.eg.db)
  library(GSVA)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
})
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE)
out<-file.path(root,"13_methodological_reanalysis_v9")
mat_out<-file.path(out,"preprocessed_strict_mapping")
score_out<-file.path(out,"module_scores_strict_mapping")
dir.create(mat_out,recursive=TRUE,showWarnings=FALSE);dir.create(score_out,recursive=TRUE,showWarnings=FALSE)
set.seed(20260802)
write_csv<-function(x,name)write.csv(x,file.path(out,name),row.names=FALSE,na="")
z<-function(x)as.numeric(scale(x));num<-function(x)suppressWarnings(as.numeric(gsub("[^0-9.+-]","",as.character(x))))
modules<-readRDS(file.path(root,"04_intermediate","gene_sets","liver_modules_combined.rds"));modules<-lapply(modules,function(x)unique(toupper(trimws(as.character(x)))))

ssgsea_all<-function(expr){sets<-lapply(modules,intersect,y=rownames(expr));sets<-sets[lengths(sets)>=3];param<-GSVA::ssgseaParam(as.matrix(expr),sets,minSize=3,normalize=TRUE);as.matrix(GSVA::gsva(param,verbose=FALSE,BPPARAM=BiocParallel::SerialParam()))}

strict_map_gpl570<-function(probes){
  valid<-intersect(probes,AnnotationDbi::keys(hgu133plus2.db::hgu133plus2.db,keytype="PROBEID"))
  a<-AnnotationDbi::select(hgu133plus2.db::hgu133plus2.db,keys=valid,keytype="PROBEID",columns="SYMBOL")
  a<-a[!is.na(a$SYMBOL)&nzchar(a$SYMBOL),];a$SYMBOL<-toupper(trimws(a$SYMBOL))
  spl<-split(a$SYMBOL,a$PROBEID);nu<-vapply(spl,function(x)length(unique(x)),integer(1));keep<-names(nu)[nu==1L]
  map<-data.frame(probe_id=keep,gene_symbol=vapply(spl[keep],function(x)unique(x)[1],character(1)),stringsAsFactors=FALSE)
  list(map=map,ambiguous=sum(nu>1L),mapped_any=length(nu),unmapped=length(probes)-length(nu))
}

strict_map_gpl4372<-function(probes){
  gpl<-readRDS(file.path(root,"02_metadata","platform_maps","GPL4372","platform_annotation.rds"))
  gpl<-gpl[match(probes,gpl$ID),c("ID","EntrezGeneID")]
  ids<-as.character(gpl$EntrezGeneID);ids[is.na(gpl$EntrezGeneID)]<-NA_character_;ids<-sub("\\.0$","",ids)
  valid<-unique(ids[!is.na(ids)&nzchar(ids)])
  a<-AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,keys=valid,keytype="ENTREZID",columns="SYMBOL")
  a<-a[!is.na(a$SYMBOL)&nzchar(a$SYMBOL),];a$SYMBOL<-toupper(trimws(a$SYMBOL))
  spl<-split(a$SYMBOL,a$ENTREZID);nu<-vapply(spl,function(x)length(unique(x)),integer(1));keep_entrez<-names(nu)[nu==1L]
  keep<-!is.na(ids)&ids%in%keep_entrez
  map<-data.frame(probe_id=as.character(gpl$ID[keep]),gene_symbol=vapply(ids[keep],function(k)unique(spl[[k]])[1],character(1)),stringsAsFactors=FALSE)
  list(map=map,ambiguous=sum(!is.na(ids)&ids%in%names(nu)[nu>1L]),mapped_any=sum(!is.na(ids)&ids%in%names(nu)),unmapped=sum(is.na(ids)|!ids%in%names(nu)))
}

impute_matrix<-function(x){storage.mode(x)<-"double";if(anyNA(x)){for(i in which(rowSums(is.na(x))>0L)){idx<-is.na(x[i,]);if(!all(idx))x[i,idx]<-median(x[i,!idx],na.rm=TRUE)}};x}
collapse_strict<-function(x,map,method=c("max_iqr","max_mean","median")){
  method<-match.arg(method);common<-intersect(rownames(x),map$probe_id);m<-map[match(common,map$probe_id),];xx<-x[common,,drop=FALSE]
  groups<-split(seq_along(common),m$gene_symbol);genes<-names(groups);ans<-matrix(NA_real_,nrow=length(groups),ncol=ncol(xx),dimnames=list(genes,colnames(xx)))
  for(j in seq_along(groups)){
    idx<-groups[[j]]
    if(length(idx)==1L)ans[j,]<-xx[idx,]
    else if(method=="max_iqr")ans[j,]<-xx[idx[which.max(apply(xx[idx,,drop=FALSE],1,IQR,na.rm=TRUE))],]
    else if(method=="max_mean")ans[j,]<-xx[idx[which.max(rowMeans(xx[idx,,drop=FALSE],na.rm=TRUE))],]
    else ans[j,]<-matrixStats::colMedians(xx[idx,,drop=FALSE],na.rm=TRUE)
  }
  vars<-apply(ans,1,var);ans[is.finite(vars)&vars>0,,drop=FALSE]
}

ds<-data.frame(dataset_id=c("MDD_GSE98793","AD_GSE48350","AD_GSE5281","AD_GSE33000"),platform=c("GPL570","GPL570","GPL570","GPL4372"),
               raw_probe_path=c("04_intermediate/expression_matrices/MDD_GSE98793/expression_matrix_raw.rds","04_intermediate/expression_matrices/AD_GSE48350/expression_matrix_raw.rds","04_intermediate/expression_matrices/AD_GSE5281/expression_matrix_raw.rds","04_intermediate/expression_matrices/AD_GSE33000/expression_matrix_raw.rds"),
               legacy_map_path=c("04_intermediate/annotation_tables/MDD_GSE98793/probe_gene_mapping_used.csv","04_intermediate/annotation_tables/AD_GSE48350/probe_gene_mapping_used.csv","04_intermediate/annotation_tables/AD_GSE5281/probe_gene_mapping_used.csv","04_intermediate/annotation_tables/AD_GSE33000/probe_gene_mapping_used.csv"),stringsAsFactors=FALSE)

pheno_path<-function(id)file.path(root,"04_intermediate","phenotype_tables",id,"pheno_harmonized.rds")
ad_map<-read.csv(file.path(root,"08_donor_aware_reanalysis","02_GSE48350_donor_region_map.csv"),check.names=FALSE)
extract<-function(fit,term){s<-coef(summary(fit));pcol<-grep("Pr\\(",colnames(s),value=TRUE)[1];c(beta=s[term,"Estimate"],SE=s[term,"Std. Error"],P=s[term,pcol])}
fit_effects<-function(id,s){ph<-readRDS(pheno_path(id));common<-intersect(colnames(s),ph$sample_id);s<-s[,common,drop=FALSE];ph<-ph[match(common,ph$sample_id),,drop=FALSE]
  do.call(rbind,lapply(rownames(s),function(m){d<-data.frame(score_z=z(s[m,]))
    if(id=="AD_GSE48350"){mp<-ad_map[match(common,ad_map$GSM),];d$diagnosis<-factor(mp$diagnosis,levels=c("Control","AD"));d$region<-factor(mp$region);d$age_z<-z(mp$age);d$sex<-factor(tolower(mp$sex));d$donor_id<-factor(mp$donor_id);fit<-lmer(score_z~diagnosis+region+age_z+sex+(1|donor_id),data=d,REML=TRUE,control=lmerControl(optimizer="bobyqa"));e<-extract(fit,"diagnosisAD");form<-"score_z ~ diagnosis + region + age_z + sex + (1 | donor_id)"
    }else if(id=="AD_GSE33000"){d$diagnosis<-factor(ph$diagnosis,levels=c("Control","AD"));d$age_z<-z(num(ph[["age:ch2"]]));d$sex<-factor(tolower(ph[["gender:ch2"]]));d<-d[!is.na(d$diagnosis),];fit<-lm(score_z~diagnosis+age_z+sex,data=d);e<-extract(fit,"diagnosisAD");form<-"score_z ~ diagnosis + age_z + sex"
    }else if(id=="AD_GSE5281"){d$diagnosis<-factor(ph$diagnosis,levels=c("Control","AD"));d$region<-factor(ph$brain_region);fit<-lm(score_z~diagnosis+region,data=d);e<-extract(fit,"diagnosisAD");form<-"descriptive score_z ~ diagnosis + region"
    }else{d$diagnosis<-factor(ph$diagnosis,levels=c("Control","MDD"));fit<-lm(score_z~diagnosis,data=d);e<-extract(fit,"diagnosisMDD");form<-"score_z ~ diagnosis"}
    data.frame(dataset_id=id,module=m,beta=e["beta"],SE=e["SE"],P=e["P"],sample_n=nrow(d),formula=form)
  }))}

annotation_rows<-list();coverage_rows<-list();result_rows<-list();summary_rows<-list()
for(i in seq_len(nrow(ds))){
  d<-ds[i,];x<-impute_matrix(readRDS(file.path(root,d$raw_probe_path)));legacy_map<-read.csv(file.path(root,d$legacy_map_path),check.names=FALSE)
  sm<-if(d$platform=="GPL570")strict_map_gpl570(rownames(x))else strict_map_gpl4372(rownames(x))
  strict_mats<-list(strict_max_iqr=collapse_strict(x,sm$map,"max_iqr"),strict_max_mean=collapse_strict(x,sm$map,"max_mean"),strict_median=collapse_strict(x,sm$map,"median"))
  legacy_mat<-readRDS(file.path(root,"13_methodological_reanalysis_v9","preprocessed_main_old_mapping",paste0(d$dataset_id,".rds")))
  mats<-c(list(legacy_first_symbol_max_mean=legacy_mat),strict_mats)
  annotation_rows[[d$dataset_id]]<-data.frame(dataset_id=d$dataset_id,platform=d$platform,initial_probe_n=nrow(x),legacy_mapped_probe_n=sum(rownames(x)%in%legacy_map$probe_id),
                                              strict_mapped_any_probe_n=sm$mapped_any,strict_ambiguous_probe_n=sm$ambiguous,strict_unmapped_probe_n=sm$unmapped,
                                              strict_retained_unambiguous_probe_n=nrow(sm$map),legacy_final_gene_n=nrow(legacy_mat),strict_IQR_final_gene_n=nrow(strict_mats$strict_max_iqr),
                                              strict_mean_final_gene_n=nrow(strict_mats$strict_max_mean),strict_median_final_gene_n=nrow(strict_mats$strict_median),
                                              ambiguity_definition="probe maps to >1 unique official gene symbol in Bioconductor annotation",stringsAsFactors=FALSE)
  scores<-list();effects<-list()
  for(method in names(mats)){
    scores[[method]]<-if(method=="legacy_first_symbol_max_mean")readRDS(file.path(root,"13_methodological_reanalysis_v9","module_scores_main_old_mapping",paste0(d$dataset_id,".rds")))else ssgsea_all(mats[[method]])
    effects[[method]]<-fit_effects(d$dataset_id,scores[[method]]);effects[[method]]$BH_FDR<-p.adjust(effects[[method]]$P,"BH")
    coverage_rows[[paste(d$dataset_id,method)]]<-data.frame(dataset_id=d$dataset_id,mapping_method=method,module=names(modules),original_module_size=lengths(modules),
                                                            covered_gene_n=vapply(modules,function(g)length(intersect(g,rownames(mats[[method]]))),integer(1)),stringsAsFactors=FALSE)
    if(method=="strict_max_iqr"){saveRDS(mats[[method]],file.path(mat_out,paste0(d$dataset_id,".rds")));saveRDS(scores[[method]],file.path(score_out,paste0(d$dataset_id,".rds")))}
  }
  base_eff<-effects[["legacy_first_symbol_max_mean"]];base_sc<-scores[["legacy_first_symbol_max_mean"]]
  for(method in setdiff(names(mats),"legacy_first_symbol_max_mean")){
    cmp_eff<-effects[[method]]
    mm0<-merge(base_eff,cmp_eff,by=c("dataset_id","module"),suffixes=c("_legacy","_comparator"))
    mm<-data.frame(dataset_id=mm0$dataset_id,module=mm0$module,comparison_mapping=method,
                   beta_legacy=mm0$beta_legacy,SE_legacy=mm0$SE_legacy,P_legacy=mm0$P_legacy,BH_FDR_legacy=mm0$BH_FDR_legacy,
                   beta_comparator=mm0$beta_comparator,SE_comparator=mm0$SE_comparator,P_comparator=mm0$P_comparator,BH_FDR_comparator=mm0$BH_FDR_comparator,
                   sample_n=mm0$sample_n_legacy,formula=mm0$formula_legacy,stringsAsFactors=FALSE)
    mm$module_score_correlation<-vapply(mm$module,function(m)cor(base_sc[m,],scores[[method]][m,],use="pairwise.complete.obs"),numeric(1))
    mm$FDR_status_change<-(mm$BH_FDR_legacy<.05)!=(mm$BH_FDR_comparator<.05)
    result_rows[[paste(d$dataset_id,method)]]<-mm
    summary_rows[[paste(d$dataset_id,method)]]<-data.frame(dataset_id=d$dataset_id,comparison_mapping=method,
      mean_score_correlation=mean(mm$module_score_correlation),min_score_correlation=min(mm$module_score_correlation),
      beta_correlation=cor(mm$beta_legacy,mm$beta_comparator),median_SE_ratio=median(mm$SE_comparator/mm$SE_legacy),FDR_status_changes=sum(mm$FDR_status_change),stringsAsFactors=FALSE)
  }
  rm(x,legacy_map,sm,strict_mats,legacy_mat,mats,scores,effects,base_eff,base_sc);invisible(gc())
}
annotation<-do.call(rbind,annotation_rows);coverage<-do.call(rbind,coverage_rows);results<-do.call(rbind,result_rows);summary<-do.call(rbind,summary_rows)
write_csv(annotation,"probe_annotation_audit.csv");write_csv(coverage,"module_gene_coverage_by_mapping.csv");write_csv(results,"probe_mapping_sensitivity_results.csv");write_csv(summary,"probe_mapping_sensitivity_summary.csv")

iqr<-results[results$comparison_mapping=="strict_max_iqr",];p<-ggplot(iqr,aes(beta_legacy,beta_comparator,color=dataset_id))+
  geom_hline(yintercept=0,color="#BDBDBD",linewidth=.3)+geom_vline(xintercept=0,color="#BDBDBD",linewidth=.3)+geom_abline(slope=1,intercept=0,linetype=2,color="#777777",linewidth=.35)+geom_point(size=1.2,alpha=.8)+facet_wrap(~dataset_id,scales="free")+
  scale_color_manual(values=c("AD_GSE33000"="#3182BD","AD_GSE48350"="#33B5A5","AD_GSE5281"="#E28E2C","MDD_GSE98793"="#D24B40"),guide="none")+
  theme_classic(base_size=7,base_family="Arial")+labs(x="Legacy first-symbol effect (SD)",y="Strict unambiguous/IQR effect (SD)",title="Probe-mapping sensitivity")
ggsave(file.path(out,"Figure_probe_mapping_beta_comparison.pdf"),p,width=183/25.4,height=120/25.4,device=grDevices::cairo_pdf,family="Arial")
ggsave(file.path(out,"Figure_probe_mapping_beta_comparison.png"),p,width=183/25.4,height=120/25.4,dpi=600,device=ragg::agg_png)

primary<-summary[summary$comparison_mapping=="strict_max_iqr"&summary$dataset_id%in%c("AD_GSE48350","AD_GSE33000"),]
stable<-all(primary$beta_correlation>=.95)&all(primary$mean_score_correlation>=.95)&all(primary$FDR_status_changes==0)
writeLines(c("# Probe mapping decision","",paste0("Decision: ",ifelse(stable,"Legacy first-symbol mapping retained as the main analysis; strict unambiguous/IQR mapping is a sensitivity analysis.","Strict unambiguous/IQR mapping is upgraded to the main analysis.")),"",
             "Pre-specified stability rule: both primary AD datasets must have beta correlation >=0.95, mean module-score correlation >=0.95, and no module-level BH-FDR transition.","",paste(capture.output(print(primary,row.names=FALSE)),collapse="\n")),file.path(out,"probe_mapping_decision.md"))
cat("Strict mapping stable for primary AD datasets:",stable,"\n");print(summary)
