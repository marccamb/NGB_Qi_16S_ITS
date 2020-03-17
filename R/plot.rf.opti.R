


plot.rf.opti <- function(foo, xlim = NULL, ylim = NULL, pch=c(22,21,24,15,16,17),
                         hue=c("#9a394e","#d68157","#f3d577","#84b368","#00876c"),
                         default.legend = F,
                         pdf.output = F, filename = NULL) {
  if (length(foo) != length(pch)) warning("The pch vector and the file vector do not have the same lenght.")

  if (pdf.output) pdf(filename, 5, 5)
  # Setting plot parameters
  par(mfrow=c(1,1), bty="l",las=1, mar=c(4,5,2,1), col.axis="gray", cex.lab=1.5)

  # Drawing an empty plot
  plot(if (is.null(xlim)) c(0,1) else xlim,
       if (is.null(ylim)) c(0,1) else ylim,
       type="n", xlab="Mean precision", ylab="Mean sensitivity")
  res_err_rate <- res_err_rate_sd <- NULL

  err <- mapply(function(files, points.type) {
    d <- readRDS(files)
    err <- mapply(function(x, col) {
      points(sensitivity_mean~precision_mean,
             pch=points.type,
             col=adjustcolor("lightgray", alpha.f = 0.3),
             data=x[-which.min(x[,"error_mean"]),])

      x.min <- x[which.min(x[,"error_mean"]),]
      segments(x.min["precision_mean"] - x.min["precision_sd"],
               x.min["sensitivity_mean"],
               x.min["precision_mean"] + x.min["precision_sd"],
               x.min["sensitivity_mean"],
               col=adjustcolor("gray", alpha.f = 0.4))

      segments(x.min["precision_mean"],
               x.min["sensitivity_mean"] - x.min["sensitivity_sd"],
               x.min["precision_mean"],
               x.min["sensitivity_mean"] + x.min["sensitivity_sd"],
               col=adjustcolor("gray", alpha.f = 0.4))
      return(list("mean"=x.min["error_mean"], "sd"=x.min["error_sd"]))
    }, d)
  }, foo, pch)

  # Points the minimum error rate for each taxo level and dataset
  mapply(function(files, points.type) {
    d <- readRDS(files)
    mapply(function(x, col){
      x.min <- x[which.min(x[,"error_mean"]),]
      points(x.min["sensitivity_mean"]~x.min["precision_mean"], cex=2,
             pch=points.type, bg="white",
             col=col)
    }, d, hue)
  },foo, pch)

  if (default.legend) {
  legend("topleft", pch=pch, col=1, bty = "n",
         legend=c("Abundant ASVs, bacteria + fungi",
                  "Abundant ASVs, bacteria",
                  "Abundant ASVs, fungi",
                  "All ASVs, bacteria + fungi",
                  "All ASVs, bacteria",
                  "All ASVs, fungi"),
         cex = 0.7)

  legend("left", pch=16, col=hue,bty = "n",
         legend=c("ASV", "Genus", "Family", "Class", "Order"),
         cex = 0.7)
  }
  if (pdf.output) dev.off()
  return(err)
}
