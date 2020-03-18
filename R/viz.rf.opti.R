#' Visualization of random forest optimization
#'
#' Visualization of the results obtained via \code{\link{rf.opti.mtry.taxo}}. Plot the mean (+/- sd) sensitivity
#'  and mean (+/- sd) precision of each models obtained for each set of parameters and/or dataset
#'
#' @param foo A vector of full names (including their path) of RDS files containing outputs of
#' the function \code{\link{rf.opti.mtry.taxo}}.
#' @param xlim A vector of lenght two giving the range of the x-axis. Dafault is c(0,1).
#' @param ylim A vector of lenght two giving the range of the y-axis. Dafault is c(0,1).
#' @param pch A vector of the same lenght as foo containing pch for plotting. Default is c(22,21,24,15,16,17).
#' @param hue A vector containing colors for each taxnonomic levels. Default is a diverging color
#'   palette of length 5.
#' @param default.legend Whether to display the default legend. Plotting your own legend is
#'   recommended. Default = F.
#' @param pdf.output Whether to save the plot in a pdf file. Default = F.
#' @param filename The filename and path where to save the pdf plot (only meaningful when pdf.output = T)
#'
#' @return Returns ....
#'
#' @export viz.rf.opti
#' @importFrom grDevices dev.off adjustcolor pdf
#' @importFrom graphics legend par plot points segments
#'

viz.rf.opti <- function(foo, xlim = c(0,1), ylim = c(0,1), pch=c(22,21,24,15,16,17),
                         hue=c("#9a394e","#d68157","#f3d577","#84b368","#00876c"),
                         default.legend = F,
                         pdf.output = F, filename = NULL) {
  if (length(foo) != length(pch)) warning("The pch vector and the file vector do not have the same lenght.")

  if (pdf.output) pdf(filename, 5, 5)
  # Setting plot parameters
  par(mfrow=c(1,1), bty="l",las=1, mar=c(4,5,2,1), col.axis="gray", cex.lab=1.5)

  # Drawing an empty plot
  plot(xlim, ylim,
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
