#' Random forest optimisation
#'
#' Runs random forest classification with several taxonomic level and mtry parameters and performs
#' k-fold or blind cross-validation.
#'
#' @param tab An abundance table containing samples in columns and OTUs/ASV in rows.
#' @param tax.table A table containing the taxonomy of each ASV/OTU.
#' @param treat A boolean vector containing the class identity of each sample, i.e. the treatment to predict.
#' This means that you should pick a class as a reference for the calculation of precision and sensitivity.
#' @param tax.lvl A character vector containing the names of the taxonomic levels to be used for asv table
#' aggregation. Default is \code{c("ASV", "genus", "family", "order", "class")}.
#' @param n.mtry The number of mtry parameters to be tested. mtry values are then calculated as
#' \code{1:n.mtry*(ncol(tab)-1)/n.mtry}. Default is 5.
#' @param cross.val The type of cross validation to perform. Possible values are "blind" or
#'  "kfold" (Default).
#' @param train.id A string that matches the name of samples tu be used for training. Only
#' meaningful for \code{cross.val = "blind"}.
#' @param n.tree The number of tree to grow for each forest. Default is 100.
#' @param cross.param The parameter needed for cross validation: the number of folds for
#' \code{cross.val = "kfold"} or the number of forests to grow for \code{cross.val = "blind"}. Default is 5.
#' @param seed The seed to set before growing each forest, and before sampling of training dataset in
#' \code{cross.val = "kfold"}. Set to NA for no seeding. Default is \code{1409}.
#' @param RDSfile A string contaning the name of the RDS file to save the results.
#' Default is NULL and results are not saved.
#'
#'@return Returns a list of dataframes corresponding to the different taxonomic levels. Each dataframe contains
#' the confusion matrix, sensitivity, precision and error rate obtained for each value of the mtry parameter.
#' Mean value and standard deviation are computed over the results of the \code{cross.param} forests grown.
#'
#' @import ranger
#' @export rf.opti.mtry.taxo
#
# 2020-02-27
# Marine C. Cambon

rf.opti.mtry.taxo <- function(tab,
                              tax.table,
                              treat,
                              n.mtry = 5,
                              tax.lvl = c("ASV", "genus", "family", "order", "class"),
                              cross.val = "kfold",
                              train.id = NA,
                              n.tree = 100,
                              cross.param = 5,
                              seed = 1409,
                              RDSfile = NULL) {

  if(!cross.val %in% c("kfold", "blind")) stop("possible methods for cross.val are \"kfold\" or \"blind\"")
  message("Ranger optimisation starting without taxonomic aggregation of the data...")
  res_tot <- list()
  for (l in tax.lvl) {
    if(l %in% c("ASV","OTU")) {
      tab_agg <- tab
    } else {
      tab_agg <- agg.table.taxo(tab, tax.lvl = l, tax.table)
    }

    ## mtry for the given taxonomic level
    if (n.mtry+1>nrow(tab_agg)) {
      n.mtry <- nrow(tab_agg)-1
      message("n.mtry is higher than the number of ", l, " in the dataset, and have been set to ", n.mtry)
    }

    res <- NULL
    for (n in 1:n.mtry) {
      if (cross.val == "kfold") tmp <- rf.kfold(tab_agg, treat,
                                                mtry = function(x) n*x/n.mtry,
                                                k.fold = cross.param,
                                                n.tree = n.tree,
                                                seed=seed)
      if (cross.val == "blind") tmp <- rf.blind(tab_agg, treat, train.id = train.id,
                                                mtry = function(x) n*x/n.mtry,
                                                n.forest = cross.param,
                                                n.tree = n.tree)
      res <- rbind(res, c(n*nrow(tab_agg)-1/n.mtry,
                          tmp[["summary"]]["mean",],
                          tmp[["summary"]]["sd",]))
    }
    colnames(res) <- c("mtry", paste(colnames(tmp[["summary"]]), "mean", sep="_"), paste(colnames(tmp[["summary"]]), "sd", sep="_"))
    res_tot[[l]] <- res

    message(l, " lvl is done\n")
  }
  if(!is.null(RDSfile)) saveRDS(res_tot, file = RDSfile)
  return(res_tot)
}
