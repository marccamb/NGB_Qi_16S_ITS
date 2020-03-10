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
res <- rf.nfold(tab, treat=z$irrigation, n_fold = 5, mtry=540, seed = 1409)
toto <-res
par(mfrow=c(2,3), mar=c(3,7,3,0), col.axis="black")
for (i in 1:length(toto$importance)) {
  barplot(sort(toto$importance[[i]], decreasing = T)[1:20],
          main=paste("n-fold", i),
          cex.names=1, cex.axis = 1.5,
          horiz = T,las=1, border=NA,
          xlab = "Gini index")
}
