path="D:\\ResearchWork\\StatisticalGenetics\\JumpSeq\\JumpSeq-project\\data\\center_three_jump\\"
aa=read.table(paste(path, "48ng_11rep_merge.sorted_GSM882244_FDR_0.0484_aa_plus_with_1_coverage.ave.perbase", sep=""), header=F)[[1]]
ac=read.table(paste(path, "48ng_11rep_merge.sorted_GSM882244_FDR_0.0484_ac_plus_with_1_coverage.ave.perbase", sep=""),header=F)[[1]]
ag=read.table(paste(path, "48ng_11rep_merge.sorted_GSM882244_FDR_0.0484_ag_plus_with_1_coverage.ave.perbase", sep=""),header=F)[[1]]
at=read.table(paste(path, "48ng_11rep_merge.sorted_GSM882244_FDR_0.0484_at_plus_with_1_coverage.ave.perbase", sep=""), header=F)[[1]]

ca=read.table(paste(path, "48ng_11rep_merge.sorted_GSM882244_FDR_0.0484_ca_plus_with_1_coverage.ave.perbase", sep=""), header=F)[[1]]
cc=read.table(paste(path, "48ng_11rep_merge.sorted_GSM882244_FDR_0.0484_cc_plus_with_1_coverage.ave.perbase", sep=""), header=F)[[1]]
cg=read.table(paste(path, "48ng_11rep_merge.sorted_GSM882244_FDR_0.0484_cg_plus_with_1_coverage.ave.perbase", sep=""), header=F)[[1]]
ct=read.table(paste(path, "48ng_11rep_merge.sorted_GSM882244_FDR_0.0484_ct_plus_with_1_coverage.ave.perbase", sep=""), header=F)[[1]]

ga=read.table(paste(path, "48ng_11rep_merge.sorted_GSM882244_FDR_0.0484_ga_plus_with_1_coverage.ave.perbase", sep=""), header=F)[[1]]
gc=read.table(paste(path, "48ng_11rep_merge.sorted_GSM882244_FDR_0.0484_gc_plus_with_1_coverage.ave.perbase", sep=""), header=F)[[1]]
gg=read.table(paste(path, "48ng_11rep_merge.sorted_GSM882244_FDR_0.0484_gg_plus_with_1_coverage.ave.perbase", sep=""),  header=F)[[1]]
gt=read.table(paste(path, "48ng_11rep_merge.sorted_GSM882244_FDR_0.0484_gt_plus_with_1_coverage.ave.perbase", sep=""), header=F)[[1]]

ta=read.table(paste(path, "48ng_11rep_merge.sorted_GSM882244_FDR_0.0484_ta_plus_with_1_coverage.ave.perbase", sep=""), header=F)[[1]]
tc=read.table(paste(path, "48ng_11rep_merge.sorted_GSM882244_FDR_0.0484_tc_plus_with_1_coverage.ave.perbase", sep=""), header=F)[[1]]
tg=read.table(paste(path, "48ng_11rep_merge.sorted_GSM882244_FDR_0.0484_tg_plus_with_1_coverage.ave.perbase", sep=""), header=F)[[1]]
tt=read.table(paste(path, "48ng_11rep_merge.sorted_GSM882244_FDR_0.0484_tt_plus_with_1_coverage.ave.perbase", sep=""), header=F)[[1]]
 
base.comp=matrix(nrow=16, ncol=19)
base.comp[1,]=(aa/sum(aa))[91:109]
base.comp[2,]=(ac/sum(ac))[91:109]
base.comp[3,]=(ag/sum(ag))[91:109]
base.comp[4,]=(at/sum(at))[91:109]
base.comp[5,]=(ca/sum(ca))[91:109]
base.comp[6,]=(cc/sum(cc))[91:109]
base.comp[7,]=(cg/sum(cg))[91:109]
base.comp[8,]=(ct/sum(ct))[91:109]
base.comp[9,]=(ga/sum(ga))[91:109]
base.comp[10,]=(gc/sum(gc))[91:109]
base.comp[11,]=(gg/sum(gg))[91:109]
base.comp[12,]=(gt/sum(gt))[91:109]
base.comp[13,]=(ta/sum(ta))[91:109]
base.comp[14,]=(tc/sum(tc))[91:109]
base.comp[15,]=(tg/sum(tg))[91:109]
base.comp[16,]=(tt/sum(tt))[91:109]

par(mfrow=c(1,2))
plot(seq(-9, 9), base.comp[1,], ylim=c(0,0.11),col=1, type="l", xlab="Distance", ylab="Reads frequency")
for (i in 2:8)
 lines(seq(-9, 9), base.comp[i,], col=i)
legend(5,0.1,  c("a*a", "a*c", "a*g", "a*t", "c*a", "c*c", "c*g", "c*t"), col=c(1,2,3,4, 5, 6, 7, 8), lty=c(1,1,1,1, 1,1,1,1))


plot(seq(-9, 9), base.comp[9,], ylim=c(0,0.11),col=1, type="l", xlab="Distance", ylab="Reads frequency")
for (i in 10:16)
 lines(seq(-9, 9), base.comp[i,], col=i)
legend(5,0.1,  c("g*a", "g*c", "g*g", "g*t", "t*a", "t*c", "t*g", "t*t"), col=c(9,10,11,12, 13, 14, 15, 16), lty=c(1,1,1,1, 1,1,1,1))
