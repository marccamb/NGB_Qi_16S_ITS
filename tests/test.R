rm(list=ls())
tmp <- readRDS("~/Documents/Postdoc_Biogeco/1_NGB/data/ORPHEE_leaf_microbiota/data_Qi_Qr_16S_filtered_formatted.RDS")
asv_table_16S <- tmp[["asv_table"]]
samples_16S <- tmp[["samples"]]
assign_16S <- tmp[["assign"]]
wp_16S <- tmp[["water_potential"]]
rm(tmp)

# Remove Qr samples from the dataset:
asv_table_16S <- asv_table_16S[,-grep("PU2", colnames(asv_table_16S))]
# Remove ASVs not present in samples:
asv_table_16S <- asv_table_16S[apply(asv_table_16S,1,sum)>0,]
# Remove assignation of ASVs not present in samples:
assign_16S <- assign_16S[match(rownames(asv_table_16S),
                               rownames(assign_16S)),]
# Remove Qr samples from the metadata
samples_16S <- samples_16S[match(colnames(asv_table_16S),
                                 rownames(samples_16S)),]

# remove endophytes samples from the dataset:
endo_16S <- subset(samples_16S, sample_type=="endo")
asv_table_endo_16S <- asv_table_16S[,match(rownames(endo_16S),
                                           dimnames(asv_table_16S)[[2]])]

samples_16S <- subset(samples_16S, sample_type=="micro")
asv_table_16S <- asv_table_16S[,match(rownames(samples_16S),
                                      dimnames(asv_table_16S)[[2]])]

## Loading ITS data -------------
tmp <- readRDS("~/Documents/Postdoc_Biogeco/1_NGB/data/ORPHEE_leaf_microbiota/data_Qi_Qr_ITS_filtered_formatted.RDS")
asv_table_ITS <- tmp[["asv_table"]]
samples_ITS <- tmp[["samples"]]
assign_ITS <- tmp[["assign"]]
wp_ITS <- tmp[["water_potential"]]
rm(tmp)

# Remove Qr samples from the dataset:
asv_table_ITS <- asv_table_ITS[,-grep("PU2", colnames(asv_table_ITS))]
# Remove ASVs not present in samples:
asv_table_ITS <- asv_table_ITS[apply(asv_table_ITS,1,sum)>0,]
# Remove assignation of ASVs not present in samples:
assign_ITS <- assign_ITS[match(rownames(asv_table_ITS),
                               rownames(assign_ITS)),]
# Remove Qr samples from the metadata
samples_ITS <- samples_ITS[match(colnames(asv_table_ITS),
                                 rownames(samples_ITS)),]

# remove endophytes samples from the dataset:
endo_ITS <- subset(samples_ITS, sample_type=="endo")
asv_table_endo_ITS <- asv_table_ITS[,match(rownames(endo_ITS),
                                           dimnames(asv_table_ITS)[[2]])]

samples_ITS <- subset(samples_ITS, sample_type=="micro")
asv_table_ITS <- asv_table_ITS[,match(rownames(samples_ITS),
                                      dimnames(asv_table_ITS)[[2]])]

## Combining 16S and ITS datasets
foo <- c(colnames(asv_table_16S), colnames(asv_table_ITS))
foo <- foo[duplicated(foo)]
tmp_16S <- asv_table_16S[,foo]
rownames(tmp_16S) <- paste(rownames(tmp_16S), "16S", sep="_")
tmp_ITS <- asv_table_ITS[,foo]
rownames(tmp_ITS) <- paste(rownames(tmp_ITS), "ITS", sep="_")
asv_table_16S_ITS <- rbind(tmp_16S, tmp_ITS)
tmp_16S <- assign_16S
rownames(tmp_16S) <- paste(rownames(tmp_16S), "16S", sep="_")
tmp_ITS <- assign_ITS
rownames(tmp_ITS) <- paste(rownames(tmp_ITS), "ITS", sep="_")
assign_16S_ITS <- rbind(tmp_16S, tmp_ITS)
samples_16S_ITS <- samples_16S[match(foo, rownames(samples_16S)),]
rm(tmp_16S, tmp_ITS)



tab <- asv_table_16S_ITS
z <- samples_16S_ITS

res <- rf.nfold(tab, treat=z$irrigation, n.fold = 5, mtry=540, seed = 1409)
toto <-res
par(mfrow=c(2,3), mar=c(3,7,3,0), col.axis="black")
for (i in 1:length(toto$importance)) {
  barplot(sort(toto$importance[[i]], decreasing = T)[1:20],
          main=paste("n-fold", i),
          cex.names=1, cex.axis = 1.5,
          horiz = T,las=1, border=NA,
          xlab = "Gini index")
}


res <- rf.blind(tab, treat = z$irrigation, train.id = "-M-", n.forest = 10)
names(res)
length(res$importance)
res$confusion
toto <-res
par(mfrow=c(2,3), mar=c(3,7,3,0), col.axis="black")
for (i in 1:length(toto$importance)) {
  barplot(sort(toto$importance[[i]], decreasing = T)[1:20],
          main=paste("n-fold", i),
          cex.names=1, cex.axis = 1.5,
          horiz = T,las=1, border=NA,
          xlab = "Gini index")
}

tab <- asv_table_16S[apply(asv_table_16S, 1, sum)>0.001*sum(asv_table_16S),]
taxo <- assign_16S[rownames(tab),]
z <- samples_16S
res_tot <- rf.opti.mtry.taxo(tab = tab,
                             tax.table = taxo,
                             n.mtry = 20,
                             treat = z$irrigation,
                             cross.val = "nfold",
                             n.tree = 500,
                             seed = 1409, rf.param = 5)
tab_agg <- agg.table.taxo(tab,tax.lvl = "family", tax.table = taxo)
dim(tab_agg)
if (n.mtry+1>nrow(tab_agg)) n.mtry <- nrow(tab_agg)-1
mtry <- 1:n.mtry*(ncol(tab_agg)-1)/n.mtry
par(mfrow=c(1,1), bty="l",las=1, mar=c(4,5,2,1), col.axis="gray", cex.lab=1.5)
plot(c(0.6,1), c(0.6,1), type="n", xlab="Mean precision", ylab="Mean sensitivity")
d <- res
for (i in 1:length(d)) {
  points(sensitivity_mean~precision_mean,
         pch=c(0,1,2,15,16,17)[i],
         col=adjustcolor("lightgray", alpha.f = 0.3),
         data=d[[i]][-which.min(d[[i]][,"error_mean"]),])
}
for (i in 1:length(d)) {
  with(data.frame(d[[i]])[which.min(d[[i]][,"error_mean"]),],
       segments(precision_mean-precision_sd, sensitivity_mean,
                precision_mean+precision_sd, sensitivity_mean,
                col=adjustcolor("gray", alpha.f = 0.4))
  )
  with(data.frame(d[[i]])[which.min(d[[i]][,"error_mean"]),],
       segments(precision_mean, sensitivity_mean-sensitivity_sd,
                precision_mean, sensitivity_mean+sensitivity_sd,
                col=adjustcolor("gray", alpha.f = 0.4))
  )
  points(sensitivity_mean~precision_mean, cex=2,
         pch=c(22,21,24,15,16,17)[i], bg="white",
         col=c(1:5)[i],
         data=data.frame(d[[i]])[which.min(d[[i]][,"error_mean"]),])
}


legend("topleft", pch=c(0,1,2,15,16,17), col=1,bty = "n",
       legend=c("Abundant ASVs, bacteria + fungi",
                "Abundant ASVs, bacteria",
                "Abundant ASVs, fungi",
                "All ASVs, bacteria + fungi",
                "All ASVs, bacteria",
                "All ASVs, fungi"),
       cex = 0.7)

legend("left", pch=16, col=hue_div_5[1:length(d)],bty = "n",
       legend=names(d),
       cex = 0.7)





res <- rf.opti.mtry.taxo(tab, treat = z$irrigation, tax.table = assign_16S_ITS,
                         cross.val = "blind",
                         train.id="-M-")

par(mfrow=c(1,1), bty="l",las=1, mar=c(4,5,2,1), col.axis="gray", cex.lab=1.5)
plot(c(0.6,1), c(0.6,1), type="n", xlab="Mean precision", ylab="Mean sensitivity")
d <- res
for (i in 1:length(d)) {
  points(sensitivity_mean~precision_mean,
         pch=c(0,1,2,15,16,17)[i],
         col=adjustcolor("lightgray", alpha.f = 0.3),
         data=d[[i]][-which.min(d[[i]][,"error_mean"]),])
}
for (i in 1:length(d)) {
  with(data.frame(d[[i]])[which.min(d[[i]][,"error_mean"]),],
       segments(precision_mean-precision_sd, sensitivity_mean,
                precision_mean+precision_sd, sensitivity_mean,
                col=adjustcolor("gray", alpha.f = 0.4))
  )
  with(data.frame(d[[i]])[which.min(d[[i]][,"error_mean"]),],
       segments(precision_mean, sensitivity_mean-sensitivity_sd,
                precision_mean, sensitivity_mean+sensitivity_sd,
                col=adjustcolor("gray", alpha.f = 0.4))
  )
  points(sensitivity_mean~precision_mean, cex=2,
         pch=c(22,21,24,15,16,17)[i], bg="white",
         col=c(1:5)[i],
         data=data.frame(d[[i]])[which.min(d[[i]][,"error_mean"]),])
}
