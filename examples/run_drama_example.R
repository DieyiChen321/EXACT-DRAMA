#!/usr/bin/env Rscript
source(file.path("R","drama.R")); i<-file.path("examples","drama","input"); o<-file.path("examples","reproduced_output","drama"); e<-file.path("examples","drama","expected_output"); dir.create(o,recursive=TRUE,showWarnings=FALSE)
rd<-function(x)read.delim(file.path(i,x),check.names=FALSE); b<-rd("beta.tsv"); s<-rd("se.tsv"); p<-rd("female_proportions.tsv"); t<-rd("true_effects.tsv"); r<-drama(b,s,cbind(p$female_proportion,1-p$female_proportion))
rmse<-function(a,b)sqrt(mean((a-b)^2)); z<-data.frame(sex=c("female","male"),original_rmse=c(rmse(b$female,t$beta.female),rmse(b$male,t$beta.male)),drama_refined_rmse=c(rmse(r$beta.female,t$beta.female),rmse(r$beta.male,t$beta.male)))
z$percent_rmse_reduction<-100*(z$original_rmse-z$drama_refined_rmse)/z$original_rmse; stopifnot(all(z$drama_refined_rmse<z$original_rmse))
write.table(r,file.path(o,"refined_effects.tsv"),sep="\t",quote=FALSE,row.names=FALSE); write.table(z,file.path(o,"accuracy_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
for(n in c("refined_effects.tsv","accuracy_summary.tsv")){a<-read.delim(file.path(e,n),check.names=FALSE);b<-read.delim(file.path(o,n),check.names=FALSE);stopifnot(isTRUE(all.equal(a,b,tolerance=1e-10,check.attributes=FALSE)))}
print(z); message("DRAMA reviewer example passed.")
