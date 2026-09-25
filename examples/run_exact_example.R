#!/usr/bin/env Rscript
if(!requireNamespace("JOBS",quietly=TRUE)) stop("Install dependencies with Rscript scripts/install_dependencies.R")
suppressPackageStartupMessages(library(MASS))
i<-file.path("examples","exact","input"); o<-file.path("examples","reproduced_output","exact"); e<-file.path("examples","exact","expected_output"); dir.create(o,recursive=TRUE,showWarnings=FALSE)
rd<-function(x) read.delim(file.path(i,x),check.names=FALSE); wr<-function(x,n) write.table(x,file.path(o,n),sep="\t",quote=FALSE,row.names=FALSE)
cb<-rd("sex_combined_beta.tsv"); cs<-rd("sex_combined_se.tsv"); fb<-rd("female_beta.tsv"); fs<-rd("female_se.tsv"); mb<-rd("male_beta.tsv"); ms<-rd("male_se.tsv"); truth<-rd("true_cell_type_effects.tsv")
w<-JOBS::jobs.nnls.weights(cb,cs); names(w)<-sub("_beta$","",names(cb)[3:ncol(cb)])
fr_all<-JOBS::jobs.eqtls(fb,fs,w,COR=FALSE); mr_all<-JOBS::jobs.eqtls(mb,ms,w,COR=FALSE); fr<-fr_all$eqtls_new; mr<-mr_all$eqtls_new; ct<-names(w)
fm<-fr[,3:ncol(fr),drop=FALSE]; mm<-mr[,3:ncol(mr),drop=FALSE]; fse<-fr_all$eqtls_se_new[,3:ncol(fr),drop=FALSE]; mse<-mr_all$eqtls_se_new[,3:ncol(mr),drop=FALSE]
fm[]<-lapply(fm,as.numeric); mm[]<-lapply(mm,as.numeric); fse[]<-lapply(fse,as.numeric); mse[]<-lapply(mse,as.numeric)
original_mean_se<-c(mean(as.matrix(fs[,3:ncol(fs)])),mean(as.matrix(ms[,3:ncol(ms)]))); refined_mean_se<-c(mean(as.matrix(fse)),mean(as.matrix(mse)))
z<-data.frame(sex=c("female","male"),original_mean_se=original_mean_se,exact_refined_mean_se=refined_mean_se)
z$percent_se_reduction<-100*(z$original_mean_se-z$exact_refined_mean_se)/z$original_mean_se; stopifnot(all(z$exact_refined_mean_se<z$original_mean_se))
wo<-data.frame(cell_type=ct,estimated_weight=as.numeric(w)); fo<-data.frame(ID=fr[[1]],fm,check.names=FALSE); mo<-data.frame(ID=mr[[1]],mm,check.names=FALSE); names(fo)[-1]<-ct; names(mo)[-1]<-ct
wr(wo,"estimated_weights.tsv"); wr(fo,"refined_female_effects.tsv"); wr(mo,"refined_male_effects.tsv"); wr(z,"precision_summary.tsv")
for(n in c("estimated_weights.tsv","refined_female_effects.tsv","refined_male_effects.tsv","precision_summary.tsv")){a<-read.delim(file.path(e,n),check.names=FALSE);b<-read.delim(file.path(o,n),check.names=FALSE);stopifnot(isTRUE(all.equal(a,b,tolerance=1e-10,check.attributes=FALSE)))}
print(z)
message("EXACT reviewer example passed.")
