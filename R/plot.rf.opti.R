

foo <- list.files(pattern = "rf_irr.*_taxo_(16S|ITS).*\\.RDS")
foo <- foo[-grep("nfold_3", foo)]

plot.rf.opti <- function(list) {

# pdf("graph_prec_sens_complete_taxo_mtry20.pdf", 5, 5)
par(mfrow=c(1,1), bty="l",las=1, mar=c(4,5,2,1), col.axis="gray", cex.lab=1.5)
plot(c(0.6,1), c(0.6,1), type="n", xlab="Mean precision", ylab="Mean sensitivity")
res_err_rate <- res_err_rate_sd <- NULL
for (j in 1:length(foo)) {
  d <- readRDS(foo[j])
  for (i in 1:length(d)) {
    points(sensitivity_mean~precision_mean,
           pch=c(0,1,2,15,16,17)[j],
           col=adjustcolor("lightgray", alpha.f = 0.3),
           data=d[[i]][-which.min(d[[i]]$err_mean),])
    res_err_rate <- c(res_err_rate, min(d[[i]]$err_mean))
    res_err_rate_sd <- c(res_err_rate_sd, min(d[[i]]$err_sd))
  }
}
for (j in 1:length(foo)) {
  d <- readRDS(foo[j])
  for (i in 1:length(d)) {
    with(d[[i]][which.min(d[[i]]$err_mean),],
         segments(precision_mean-precision_sd, sensitivity_mean,
                  precision_mean+precision_sd, sensitivity_mean,
                  col=adjustcolor("gray", alpha.f = 0.4))
    )
    with(d[[i]][which.min(d[[i]]$err_mean),],
         segments(precision_mean, sensitivity_mean-sensitivity_sd,
                  precision_mean, sensitivity_mean+sensitivity_sd,
                  col=adjustcolor("gray", alpha.f = 0.4))
    )
    points(sensitivity_mean~precision_mean, cex=2,
           pch=c(22,21,24,15,16,17)[j], bg="white",
           col=hue_div_5[i],
           data=d[[i]][which.min(d[[i]]$err_mean),])
  }
}

# legend("topleft", pch=c(0,1,2,15,16,17), col=1,bty = "n",
#        legend=gsub("rf_irr_(.*)\\.RDS","\\1",foo),
#        cex = 0.7)

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
}
