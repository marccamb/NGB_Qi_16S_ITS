# Funtion rf.mtry.taxo()
# runs random forest classification with several
#  taxonomic level and mtry parameters and returns
#  error rate obtained by n-foldinf cross-validation
#
# 2020-02-27
# Marine C. Cambon

# # For debugging:
# l <- "ASV"
# i <- 1
# n <- 1
# n_fold <- 5
# n.mtry <- 5
# seed <- 1409
# tax.table <- taxo
# variable <- "irrigation"
# treat <- z$irrigation
# train.id <- "-M-"

rf.opti.mtry.taxo <- function(tab, tax.table, treat,
                              n.mtry = 5,
                              cross.val = "nfold",
                              train.id = NA,
                              n.tree = 100,
                              rf.param = 5,
                              seed = 1409) {


  res_tot <- list()
  for (l in c("ASV", "genus", "family", "order", "class")) {
    if(l=="ASV") {
      tab_agg <- tab
    } else {
      tab_agg <- agg.table.taxo(tab, tax.lvl = l, taxo)
    }

    ## mtry for the given taxonomic level
    if (n.mtry+1>ncol(tab_agg)) n.mtry <- ncol(tab_agg)-1
    mtry <- 1:n.mtry*(ncol(tab_agg)-1)/n.mtry

    res <- NULL
    for (n in 1:n.mtry) {
      if (cross.val == "nfold") tmp <- rf.nfold(tab_agg, treat,
                                                mtry = function(x) n*x/n.mtry,
                                                n.fold = rf.param,
                                                n.tree = n.tree,
                                                seed=seed)
      res <- rbind(res, c(mtry[n],
                          tmp[["summary"]]["mean",],
                          tmp[["summary"]]["sd",]))
    }
    colnames(res) <- c("mtry", paste(colnames(tmp[["summary"]]), "mean", sep="_"), paste(colnames(tmp[["summary"]]), "sd", sep="_"))
    res_tot[[l]] <- res

    message(l, " lvl is done\n")
  }
  return(res_tot)
}

