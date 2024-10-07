
#univariate MR
PATH<-""
setwd(PATH)

library(ggplot2)
library(plyr)
library(data.table)
library(devtools)
library(MendelianRandomization)
library(TwoSampleMR)

EXP1<-data.table::fread("prot-a-N.clump.txt")
EXP1<-data.frame(EXP1)
head(EXP1)

OUT1<-extract_outcome_data(snps=EXP1_IV$SNP,outcomes="finn-b-I9_HYPTENSPUL",proxies=T,maf_threshold = 0.01,access_token = NULL)
OUT1<-OUT1[!duplicated(OUT1$SNP),]
OUT1$id.outcome<-"PAH"
OUT1$outcome<-"PAH"
head(OUT1)
data_h<-harmonise_data(exposure_dat=EXP1_IV,outcome_dat=OUT1,action=2)
mr_outcome<-mr(data_h)
#(univariate reverse MR)
#PATH<-""
#setwd(PATH)
#EXP<-extract_instruments(outcomes = "finn-b-I9_HYPTENSPUL", p1 = 5e-06,clump = T,r2=0.001,kb=10000, p2 = 5e-08)
#OUT1<-data.table::fread("PROT-N.gz")
#OUT<-data.frame(OUT1)
#OUT$id.outcome<-"CTSs"
#OUT$outcome<-"CTSs"
#total<-merge(OUT,EXP,by.x="SNP",by.y="SNP",all = F)
#total<-subset(total,pval.outcome>5e-08)
#total<-total[!duplicated(total$SNP),]
#EXP<-total[,c("SNP","effect_allele.exposure","other_allele.exposure", "eaf.exposure", "beta.exposure","se.exposure", "pval.exposure","id.exposure","exposure","samplesize.exposure")]
#OUT<-total[,c("SNP","effect_allele.outcome","other_allele.outcome", "eaf.outcome", "beta.outcome","se.outcome","pval.outcome","id.outcome","outcome","samplesize.outcome")]
#data_h<-harmonise_data(exposure_dat=EXP,outcome_dat=OUT,action=2)
#mr_outcome<-mr(data_h)
#mr_outcome
mr_outcome
H<-mr_heterogeneity(data_h)
H
ple <- mr_pleiotropy_test(data_h)
ple
library(MRPRESSO)
presso<-run_mr_presso(data_h, NbDistribution = 1000, SignifThreshold = 0.05)
presso


library(MendelianRandomization)
library(TwoSampleMR)
#multivariate MR
setwd("")
df=read.table("cathepsin_id.txt",header = T,sep = "\t")
id=as.vector(df$id)
expo_id=c(id)
expo_rt<- mv_extract_exposures(expo_id,pval_threshold = 5e-06)
outcome_dat <- read_outcome_data(
  snps = expo_rt$SNP,
  filename = "finngen_R9_I9_HYPTENSPUL.txt",
  sep = "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value")

mvhr_rt <- mv_harmonise_data(expo_rt, outcome_dat)

res <- mv_multiple(mvhr_rt)
res$result$or <- exp(res$result$b)
res$result$lci95 <- res$result$b-1.96*res$result$se
res$result$upci95 <- res$result$b+1.96*res$result$se
res$result$or_lci95 <- exp(res$result$lci95)
res$result$or_upci95 <- exp(res$result$upci95)

write.table(res$result,"OR.txt",row.names = F,quote = F,sep = "\t")

mvmrdata <- mr_mvinput(bx = mvhr_rt$exposure_beta, bxse = mvhr_rt$exposure_se, by = mvhr_rt$outcome_beta, byse = mvhr_rt$outcome_se, correlation =matrix())
res_ivw<-mr_mvivw(mvmrdata)
res_ivw
res_ergger<-mr_mvegger(mvmrdata)
res_ergger
