########################################
# Script: 06_recompute_primary_strict.R
# Purpose: Recompute locked strict-mapping dataset effects, donor-aware comparisons, CR2 sensitivity, REML, and Hartung-Knapp summaries.
# Input: Strict-mapping module scores, phenotype data, and donor map.
# Output: Dataset-specific effects, model comparisons, meta-analysis tables, and figures.
# Software: R
# Version: 4.5.2
# Random seed: Not applicable (deterministic)
# Author: Study authors
########################################
options(stringsAsFactors=FALSE,width=180,contrasts=c("contr.treatment","contr.poly"))
suppressPackageStartupMessages({library(lme4);library(lmerTest);library(metafor);library(ggplot2)})
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE);out<-file.path(root,"13_methodological_reanalysis_v9")
set.seed(20260802)
write_csv<-function(x,name)write.csv(x,file.path(out,name),row.names=FALSE,na="")
z<-function(x)as.numeric(scale(x));num<-function(x)suppressWarnings(as.numeric(gsub("[^0-9.+-]","",as.character(x))))
extract<-function(fit,term){s<-coef(summary(fit));if(!term%in%rownames(s))return(c(beta=NA,SE=NA,df=NA,CI_low=NA,CI_high=NA,P=NA));pcol<-grep("Pr\\(",colnames(s),value=TRUE)[1];df<-if("df"%in%colnames(s))s[term,"df"]else Inf;crit<-if(is.finite(df))qt(.975,df)else qnorm(.975);c(beta=s[term,"Estimate"],SE=s[term,"Std. Error"],df=df,CI_low=s[term,"Estimate"]-crit*s[term,"Std. Error"],CI_high=s[term,"Estimate"]+crit*s[term,"Std. Error"],P=s[term,pcol])}
status<-function(x)ifelse(is.finite(x)&x<.05,"FDR<0.05","not_FDR<0.05")
score_dir<-file.path(out,"module_scores_strict_mapping")
module_names<-rownames(readRDS(file.path(score_dir,"AD_GSE48350.rds")))

# GSE48350 exact same 253 samples, strict mapping main.
s483<-readRDS(file.path(score_dir,"AD_GSE48350.rds"));mp<-read.csv(file.path(root,"08_donor_aware_reanalysis","02_GSE48350_donor_region_map.csv"),check.names=FALSE);mp<-mp[match(colnames(s483),mp$GSM),]
d0<-data.frame(diagnosis=factor(mp$diagnosis,levels=c("Control","AD")),region=factor(mp$region),age_z=z(mp$age),sex=factor(tolower(mp$sex)),donor_id=factor(mp$donor_id))
if(nrow(d0)!=253||anyNA(d0))stop("GSE48350 strict main mapping failed sample lock")
r483<-do.call(rbind,lapply(module_names,function(m){d<-d0;d$score_z<-z(s483[m,]);nv<-lm(score_z~diagnosis+region+age_z+sex,data=d);da<-lmer(score_z~diagnosis+region+age_z+sex+(1|donor_id),data=d,REML=TRUE,control=lmerControl(optimizer="bobyqa"));en<-extract(nv,"diagnosisAD");ed<-extract(da,"diagnosisAD")
  cr<-tryCatch({ct<-clubSandwich::coef_test(nv,vcov="CR2",cluster=d$donor_id,test="Satterthwaite");rr<-ct[rownames(ct)=="diagnosisAD",];c(beta=rr$beta,SE=rr$SE,df=rr$df_Satt,P=rr$p_Satt)},error=function(e)c(beta=NA,SE=NA,df=NA,P=NA));vc<-as.data.frame(VarCorr(da))
  data.frame(module=m,naive_beta=en["beta"],naive_SE=en["SE"],naive_df=en["df"],naive_CI_low=en["CI_low"],naive_CI_high=en["CI_high"],naive_P=en["P"],donoraware_beta=ed["beta"],donoraware_SE=ed["SE"],donoraware_df=ed["df"],donoraware_CI_low=ed["CI_low"],donoraware_CI_high=ed["CI_high"],donoraware_P=ed["P"],cluster_robust_beta=cr["beta"],cluster_robust_SE=cr["SE"],cluster_robust_df=cr["df"],cluster_robust_P=cr["P"],sample_n=nrow(d),donor_n=nlevels(d$donor_id),random_intercept_variance=vc$vcov[vc$grp=="donor_id"][1],residual_variance=vc$vcov[vc$grp=="Residual"][1],singular_fit=isSingular(da,tol=1e-5),naive_formula="score_z ~ diagnosis + region + age_z + sex",donoraware_formula="score_z ~ diagnosis + region + age_z + sex + (1 | donor_id)",cluster_robust_formula="same fixed-effects OLS with donor-clustered CR2 SE",mapping="strict one-to-one official symbol; highest-IQR probe")
}))
r483$cluster_robust_CI_low<-r483$cluster_robust_beta-qt(.975,r483$cluster_robust_df)*r483$cluster_robust_SE
r483$cluster_robust_CI_high<-r483$cluster_robust_beta+qt(.975,r483$cluster_robust_df)*r483$cluster_robust_SE
r483$naive_BH_FDR<-p.adjust(r483$naive_P,"BH");r483$donoraware_BH_FDR<-p.adjust(r483$donoraware_P,"BH");r483$cluster_robust_BH_FDR<-p.adjust(r483$cluster_robust_P,"BH");r483$beta_difference_donoraware_minus_naive<-r483$donoraware_beta-r483$naive_beta;r483$SE_ratio_donoraware_over_naive<-r483$donoraware_SE/r483$naive_SE;r483$significance_transition<-paste(status(r483$naive_BH_FDR),"to",status(r483$donoraware_BH_FDR))
write_csv(r483,"GSE48350_naive_vs_donoraware_full.csv")

# GSE33000 strict mapping main.
s330<-readRDS(file.path(score_dir,"AD_GSE33000.rds"));ph330<-readRDS(file.path(root,"04_intermediate","phenotype_tables","AD_GSE33000","pheno_harmonized.rds"));ph330<-ph330[match(colnames(s330),ph330$sample_id),];keep<-ph330$diagnosis%in%c("AD","Control");s330<-s330[,keep];ph330<-ph330[keep,]
d330<-data.frame(diagnosis=factor(ph330$diagnosis,levels=c("Control","AD")),age_z=z(num(ph330[["age:ch2"]])),sex=factor(tolower(ph330[["gender:ch2"]])))
r330<-do.call(rbind,lapply(module_names,function(m){d<-d330;d$score_z<-z(s330[m,]);fit<-lm(score_z~diagnosis+age_z+sex,data=d);e<-extract(fit,"diagnosisAD");data.frame(dataset="GSE33000",module=m,beta=e["beta"],SE=e["SE"],df=e["df"],CI_low=e["CI_low"],CI_high=e["CI_high"],P=e["P"],sample_n=nrow(d),donor_n=NA,formula="score_z ~ diagnosis + age_z + sex",mapping="strict one-to-one official symbol; highest-IQR probe")}));r330$BH_FDR<-p.adjust(r330$P,"BH")
r483_primary<-data.frame(dataset="GSE48350",module=r483$module,beta=r483$donoraware_beta,SE=r483$donoraware_SE,df=r483$donoraware_df,CI_low=r483$donoraware_CI_low,CI_high=r483$donoraware_CI_high,P=r483$donoraware_P,sample_n=r483$sample_n,donor_n=r483$donor_n,formula=r483$donoraware_formula,mapping=r483$mapping,BH_FDR=r483$donoraware_BH_FDR)
ad_effects<-rbind(r483_primary,r330);write_csv(ad_effects,"AD_dataset_specific_effects_strict_main.csv")

# Conventional REML and Hartung-Knapp are kept separate.
meta_rows<-list();hk_rows<-list()
for(m in module_names){d<-ad_effects[ad_effects$module==m,];fit<-metafor::rma(yi=d$beta,sei=d$SE,method="REML");hk<-metafor::rma(yi=d$beta,sei=d$SE,method="REML",test="knha")
  meta_rows[[m]]<-data.frame(module=m,k=2,effect=as.numeric(fit$b),SE=fit$se,CI_low=fit$ci.lb,CI_high=fit$ci.ub,P=fit$pval,tau2=fit$tau2,I2=fit$I2,Q=fit$QE,method="REML normal-theory")
  hk_rows[[m]]<-data.frame(module=m,k=2,effect=as.numeric(hk$b),SE=hk$se,CI_low=hk$ci.lb,CI_high=hk$ci.ub,P=hk$pval,tau2=hk$tau2,I2=hk$I2,Q=hk$QE,method="REML with Hartung-Knapp",df=hk$ddf)
}
reml<-do.call(rbind,meta_rows);hk<-do.call(rbind,hk_rows);reml$BH_FDR<-p.adjust(reml$P,"BH");hk$BH_FDR<-p.adjust(hk$P,"BH");write_csv(reml,"AD_two_dataset_REML_strict_main.csv");write_csv(hk,"AD_two_dataset_Hartung_Knapp_strict_main.csv")

# MDD blood strict mapping, retaining the project's Welch/Cohen-d specification.
s987<-readRDS(file.path(score_dir,"MDD_GSE98793.rds"));ph987<-readRDS(file.path(root,"04_intermediate","phenotype_tables","MDD_GSE98793","pheno_harmonized.rds"));ph987<-ph987[match(colnames(s987),ph987$sample_id),]
mdd_blood<-do.call(rbind,lapply(module_names,function(m){case<-as.numeric(s987[m,ph987$diagnosis=="MDD"]);ctl<-as.numeric(s987[m,ph987$diagnosis=="Control"]);n1<-length(case);n0<-length(ctl);sp<-sqrt(((n1-1)*var(case)+(n0-1)*var(ctl))/(n1+n0-2));d<-(mean(case)-mean(ctl))/sp;se<-sqrt(1/n1+1/n0+d^2/(2*(n1+n0)));tt<-t.test(case,ctl);data.frame(module=m,effect=d,SE=se,CI_low=d-1.96*se,CI_high=d+1.96*se,P=tt$p.value,MDD_n=n1,control_n=n0,sample_n=n1+n0,donor_n=n1+n0,formula="Welch two-sample t-test; Cohen d = MDD minus Control",mapping="strict one-to-one official symbol; highest-IQR probe")}));mdd_blood$BH_FDR<-p.adjust(mdd_blood$P,"BH");write_csv(mdd_blood,"MDD_GSE98793_strict_main_results.csv")

# Updated same-data figures.
pal<-c("FDR<0.05 to not_FDR<0.05"="#D24B40","FDR<0.05 to FDR<0.05"="#3182BD","not_FDR<0.05 to FDR<0.05"="#33B5A5","not_FDR<0.05 to not_FDR<0.05"="#767676")
theme_set(theme_classic(base_size=7,base_family="Arial"))
p1<-ggplot(r483,aes(naive_beta,donoraware_beta,color=significance_transition))+geom_hline(yintercept=0,color="#BDBDBD",linewidth=.3)+geom_vline(xintercept=0,color="#BDBDBD",linewidth=.3)+geom_abline(slope=1,intercept=0,linetype=2,color="#777777",linewidth=.35)+geom_point(size=1.7)+scale_color_manual(values=pal,drop=FALSE)+coord_equal()+labs(x="Naive diagnosis coefficient (SD)",y="Donor-aware diagnosis coefficient (SD)",color="BH-FDR transition",title="GSE48350 same-sample estimates (strict mapping)")
p2<-ggplot(r483,aes(reorder(module,SE_ratio_donoraware_over_naive),SE_ratio_donoraware_over_naive,fill=significance_transition))+geom_hline(yintercept=1,linetype=2,color="#777777",linewidth=.35)+geom_col(width=.72)+coord_flip()+scale_fill_manual(values=pal,drop=FALSE)+labs(x=NULL,y="Donor-aware SE / naive SE",fill="BH-FDR transition",title="Uncertainty after donor adjustment")
tc<-as.data.frame(table(r483$significance_transition));names(tc)<-c("transition","modules");p3<-ggplot(tc,aes(reorder(transition,modules),modules,fill=transition))+geom_col(width=.65)+coord_flip()+scale_fill_manual(values=pal,guide="none")+labs(x=NULL,y="Modules (n)",title="BH-FDR transitions across 24 modules")
save<-function(p,stem,w,h){ggsave(file.path(out,paste0(stem,".pdf")),p,width=w/25.4,height=h/25.4,device=grDevices::cairo_pdf,family="Arial");ggsave(file.path(out,paste0(stem,".png")),p,width=w/25.4,height=h/25.4,dpi=600,device=ragg::agg_png)}
save(p1,"Figure_GSE48350_beta_comparison",89,80);save(p2,"Figure_GSE48350_SE_ratio",120,115);save(p3,"Figure_GSE48350_significance_transition",89,70)
cat("GSE48350 strict naive/donor-aware FDR:",sum(r483$naive_BH_FDR<.05),sum(r483$donoraware_BH_FDR<.05)," beta cor=",cor(r483$naive_beta,r483$donoraware_beta)," median SE ratio=",median(r483$SE_ratio_donoraware_over_naive),"\n")
cat("GSE33000 strict FDR:",sum(r330$BH_FDR<.05)," REML FDR:",sum(reml$BH_FDR<.05)," HK FDR:",sum(hk$BH_FDR<.05)," MDD blood FDR:",sum(mdd_blood$BH_FDR<.05),"\n")
