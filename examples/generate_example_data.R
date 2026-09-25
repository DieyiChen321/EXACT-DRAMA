#!/usr/bin/env Rscript
set.seed(20260925)
root <- "examples"
cell_types <- c("B_IN","B_Mem","Plasma","CD4_ET","CD4_NC","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","DC1","Mono_C","Mono_NC","NK","NK_R")
write_tsv <- function(x, path) { dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE); write.table(x,path,sep="\t",quote=FALSE,row.names=FALSE) }

# EXACT: precise bulk effects are weighted combinations of 14 noisy cell-type effects.
n <- 240L; ids <- sprintf("ENSGEXAMPLE%03d-rs%06d",1:n,1:n)
w <- c(.13,.09,.04,.08,.13,.05,.07,.06,.04,.03,.08,.07,.07,.06)
base <- matrix(rnorm(n*14,sd=.18),ncol=14); shift <- matrix(rnorm(n*14,sd=.035),ncol=14)
truth_f <- base+shift; truth_m <- base-shift; truth_c <- (truth_f+truth_m)/2
sc_se <- .13; bulk_se <- .025
observe <- function(x) x+matrix(rnorm(length(x),sd=sc_se),ncol=14)
obs_c <- truth_c+matrix(rnorm(length(truth_c),sd=.01),ncol=14); obs_f <- observe(truth_f); obs_m <- observe(truth_m)
bulk <- function(x) as.numeric(x%*%w+rnorm(n,sd=bulk_se))
make_beta <- function(b,s) { x<-data.frame(ID=ids,beta.est=b,s,check.names=FALSE); names(x)[3:ncol(x)]<-paste0(cell_types,"_beta"); x }
make_se <- function() { x<-data.frame(ID=ids,beta.est=rep(bulk_se,n),matrix(sc_se,n,14),check.names=FALSE); names(x)[3:ncol(x)]<-paste0(cell_types,"_beta"); x }
d <- file.path(root,"exact","input")
write_tsv(make_beta(bulk(truth_c),obs_c),file.path(d,"sex_combined_beta.tsv")); write_tsv(make_se(),file.path(d,"sex_combined_se.tsv"))
write_tsv(make_beta(bulk(truth_f),obs_f),file.path(d,"female_beta.tsv")); write_tsv(make_se(),file.path(d,"female_se.tsv"))
write_tsv(make_beta(bulk(truth_m),obs_m),file.path(d,"male_beta.tsv")); write_tsv(make_se(),file.path(d,"male_se.tsv"))
truth<-data.frame(ID=rep(ids,2),sex=rep(c("female","male"),each=n),rbind(truth_f,truth_m),check.names=FALSE); names(truth)[3:ncol(truth)]<-cell_types
write_tsv(truth,file.path(d,"true_cell_type_effects.tsv")); write_tsv(data.frame(cell_type=cell_types,true_weight=w),file.path(d,"true_weights.tsv"))

# DRAMA: noisy sex-only GWAS plus three precise sex-combined studies.
n <- 300L; snp<-sprintf("1_%d_A_G_b38",100000L+1:n); tf<-rnorm(n,sd=.16); tm<-tf+rnorm(n,sd=.07)
sex_se <- .20; cs<-c(.055,.060,.050); p<-c(1,0,.62,.48,.35)
b<-data.frame(snp=snp,female=tf+rnorm(n,sd=sex_se),male=tm+rnorm(n,sd=sex_se),
 combined_1=p[3]*tf+(1-p[3])*tm+rnorm(n,sd=cs[1]),combined_2=p[4]*tf+(1-p[4])*tm+rnorm(n,sd=cs[2]),combined_3=p[5]*tf+(1-p[5])*tm+rnorm(n,sd=cs[3]))
s<-data.frame(snp=snp,female=sex_se,male=sex_se,combined_1=cs[1],combined_2=cs[2],combined_3=cs[3])
d<-file.path(root,"drama","input"); write_tsv(b,file.path(d,"beta.tsv")); write_tsv(s,file.path(d,"se.tsv"))
write_tsv(data.frame(study=names(b)[-1],female_proportion=p),file.path(d,"female_proportions.tsv")); write_tsv(data.frame(snp=snp,beta.female=tf,beta.male=tm),file.path(d,"true_effects.tsv"))
message("Generated deterministic example inputs.")
