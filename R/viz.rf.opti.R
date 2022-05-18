#' Visualization of random forest optimization
#'
#' Visualization of the results obtained via \code{\link{rf.opti.mtry.taxo}}.
#' Plot the mean (+/- sd) sensitivity and mean (+/- sd) precision of each models
#'  obtained for each set of parameters and/or dataset
#'
#' @param foo The output of \code{\link{rf.opti.mtry.taxo}}, or alternatively,
#' a vector of full names (including their path, usually obtained with
#' \code{list.files(pattern=".*RDS", full.names=T)}) of RDS files
#' containing outputs of the function \code{\link{rf.opti.mtry.taxo}}.
#' @param plot Whether to plot the graph. Default is TRUE. If FALSE, only returns
#'  the features of the best prediction (see 'Value').
#' @param all_mtry Whether to plot the results from all mtry values. Default is
#' TRUE. If FALSE, only the results for the best mtry parameters will be plotted.
#' @param display_tax An optionnal vector of the names of the taxonomic levels
#' to plot. Default = NULL plot all levels.
#' @param xlim A vector of length two giving the range of the x-axis.
#' Default is c(0,1).
#' @param ylim A vector of length two giving the range of the y-axis.
#' Default is c(0,1).
#' @param pch A vector of the same length as foo containing pch for plotting.
#' Default is c(22,21,24,15,16,17).
#' @param hue A vector containing colors for each taxonomic levels.
#' Default is a diverging color palette of length 5.
#' @param axis.labels A a vector of length 2 containing x-axis and y-axis labels.
#' @param default.legend Whether to display the default legend. Default = T.
#' @param pdf.output Whether to save the plot in a pdf file. Default = F.
#' @param filename The file name and path where to save the pdf plot (only
#' meaningful when pdf.output = T).
#' @param ... Other graphical parameters to pass to plot().
#'
#' @return Returns the features of the best prediction (i.e. giving the lowest
#' mean error rate), including the
#' corresponding file name and taxonomic level.
#'
#' @export viz.rf.opti
#' @importFrom grDevices dev.off adjustcolor pdf
#' @importFrom graphics legend par plot points segments
#' @import purrr

viz.rf.opti <- function(foo, plot = TRUE,
                        xlim = c(0,1), ylim = c(0,1), pch=c(22,21,24,15,16,17),
                        hue=c("#9a394e","#d68157","#f3d577","#84b368","#00876c"),
                        axis.labels = c("Mean precision","Mean sensitivity"),
                        default.legend = TRUE, all_mtry=FALSE, display_tax=NULL,
                        pdf.output = FALSE, filename = NULL, ...) {
  # Check if foo contains RDS file names
  RDSfiles = F
  if (class(foo)!="list") {
    if (length(grep("\\.RDS", foo)) != length(foo)) {
    stop("foo should contain names of RDS files. Check the RDS extension of your files.")
    } else {
      RDSfiles = T
    }
  }
  # Opening a PDF file if needed
  if (plot & pdf.output) pdf(filename, 5, 5)
  # Setting plot parameters
  if (plot) {
    if(is.null(display_tax)) {
      tax.lvl <- names(readRDS(foo[1]))
    } else {
      tax.lvl <- display_tax
    }
    if (length(tax.lvl) != length(hue)) warning("The hue vector and the taxonomic levels do not have the same length.")

    # Drawing an empty plot
    par(bty="l",las=1, mar=c(4,5,2,1), col.axis="gray30", cex.lab=1.5)
    plot(xlim, ylim,
         type="n", xlab=axis.labels[1], ylab=axis.labels[2], ...)
  }

  # plot with multiple RDS files
  if (RDSfiles) {
    if (plot & length(foo) != length(pch)) warning("The pch vector and the file vector do not have the same length.")
    res_err <- purrr::map2(foo, pch, function(files, points.type) {
      d <- readRDS(files)
      # Only keeping geven tax lvl
      if (!is.null(display_tax)) d <- d[names(d) %in% display_tax]

      err <- purrr::map(d, function(x) {
        if(nrow(x)<2) {
          x.min <- x
        } else {
          x.min <- x[which.min(x[,"error_mean"]),]
        }
        if(plot) {
          if(all_mtry) {
            points(x[-which.min(x[,"error_mean"]),"sensitivity_mean"]~
                     x[-which.min(x[,"error_mean"]),"precision_mean"],
                   pch=points.type,
                   col=adjustcolor("lightgray", alpha.f = 0.3))
          }

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
        }
        return(x.min)
      })
      err <- do.call(rbind, err)
      rownames(err) <- names(d)
      return(err)
    })
    names(res_err) <- lapply(strsplit(foo, "/"), function(x) x[length(x)])

    # Points the minimum error rate for each taxo level and dataset
    if(plot) {
      purrr::walk2(res_err, pch, function(x, points.type) {
          points(x[,"sensitivity_mean"]~x[,"precision_mean"], cex=2,
                 pch=points.type, bg="white",
                 col=hue)
      })
    }
    best <- purrr::map(res_err, function(x) x[which.min(x[,"error_mean"]),])
    best <- do.call(rbind, best)
    res_err[["best"]] <- c("file"=rownames(best)[which.min(best[,"error_mean"])],
                           best[which.min(best[,"error_mean"]),])
  } else {
# Plot with only one object from rf.opti.mtry.taxo
    if (length(pch)>1) {
      warning("There is only one object to plot. Only the first value of pch will be used.")
    }
    if(is.null(display_tax)) {
      tax.lvl <- names(readRDS(foo[1]))
    } else {
      tax.lvl <- display_tax
    }
    if (length(tax.lvl) != length(hue)) warning("The hue vector and the taxonomic levels do not have the same length.")

    err <- mapply(function(x, col) {
      x <- x[apply(x, 1, function(r) !any(!is.finite(r))),]
      x.min <- x[which.min(x[,"error_mean"]),]
      if(plot) {
        points(sensitivity_mean~precision_mean,
               pch=pch[1],
               col=adjustcolor("lightgray", alpha.f = 0.3),
               data=x[-which.min(x[,"error_mean"]),])
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
      }
      return(list("mean"=x.min["error_mean"], "sd"=x.min["error_sd"]))
    }, foo)

    # Points the minimum error rate for each taxo level
    if(plot) {
      mapply(function(x, col){
        x <- x[apply(x, 1, function(r) !any(!is.finite(r))),]
        x.min <- x[which.min(x[,"error_mean"]),]
        points(x.min["sensitivity_mean"]~x.min["precision_mean"], cex=2,
               pch=pch[1], bg="white",
               col=col)
      }, foo, hue)
    }
    res_err <- err[,which.min(err["mean",])]
    names(res_err) <- paste(names(res_err), colnames(err)[which.min(err["mean",])])
  }

  if (plot & default.legend) {
    legend("topleft", pch=pch[1:length(foo)], col=1, bty = "n",
           legend=gsub("^.*/(.*)\\.RDS$", "\\1", foo),
           cex = 0.7)

    legend("left", pch=16, col=hue[1:length(tax.lvl)], bty = "n",
           legend=tax.lvl,
           cex = 0.7)
  }
  if (plot & pdf.output) dev.off()
  return(res_err)
}
