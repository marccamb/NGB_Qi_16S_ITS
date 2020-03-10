#' k-fold cross validation for random forest
#'
#' Grows \code{n} random forests for classification and \code{k}-fold cross validation (\code{k=n})
#'
#' @param tab An abundance table containing samples in columns and OTUs/ASV in rows.
#' @param treat A vector containing the class identity of each sample.
#' @param n.fold A number of fold to be applied for n-fold cross-valisation.
#' @param mtry The mtry parameter to be passed to the \code{ranger} function. See \code{ranger} documentation for details.
#' @param n.tree The number of tree to grow. The default is \code{500}.
#' @param seed A number to set the seed before before growing each forest. The default is \code{NA}.
#'
#'
#' @return A list...
#'
#' @import ranger

# 2020-02-27
# Marine C. Cambon


rf.nfold <- function(tab, treat,
                     n.fold = 5,
                     mtry = NULL,
                     n.tree = 500,
                     seed = NA) {
  # Preparing training IDs and dataframe
  if (!is.na(seed)) set.seed(seed)
  train.idx <- sample(rep(1:n.fold, 1/n.fold * ncol(tab)), replace = F)
  tab_agg <- data.frame("treat" = treat, t(tab))

  res <- NULL
  importance <- list()
  err_mean <- err_sd <- NULL
  TP_mean <- TN_mean <- FP_mean <- FN_mean <- NULL
  TP_sd <- TN_sd <- FP_sd <- FN_sd <- NULL
  sensitivity_mean <- precision_mean <- NULL
  sensitivity_sd <- precision_sd <- NULL

  for (i in 1:n.fold) {
    # Split training and test datasets
    train <- tab_agg[train.idx != i, ]
    test <- tab_agg[train.idx == i, ]

    # Grow the forest and make predictions
    if (!is.na(seed)) set.seed(seed)
    rg <- ranger::ranger(treat ~ ., data = train,
                 num.trees = n.tree,
                 mtry = mtry,
                 importance = "impurity")
    pred <- stats::predict(rg, data = test)

    # Store the variables
    tmp <- data.frame(table(pred$predictions, test$treat))
    TN <- tmp[tmp$Var1=="irr" & tmp$Var2=="irr","Freq"]
    TP <- tmp[tmp$Var1=="non-irr" & tmp$Var2=="non-irr","Freq"]
    FN <- tmp[tmp$Var1=="non-irr" & tmp$Var2=="irr","Freq"]
    FP <- tmp[tmp$Var1=="irr" & tmp$Var2=="non-irr","Freq"]
    error <- (FP+FN)/(FP+FN+TP+TN)
    sensitivity <- TP/(TP+FN)
    precision <- TP/(TP+FP)
    res <- rbind(res,c(TN,TP,FN,FP,error,sensitivity,precision))
    importance[[i]] <- ranger::importance(rg)
  }
  colnames(res) <- c("TN","TP","FN","FP","error","sensitivity","precision")
  rownames(res) <- paste("nfold_", 1:n.fold,sep="")
  names(importance) <- paste("nfold_", 1:n.fold,sep="")

  summary <- rbind(apply(res,2,mean),apply(res,2,sd))
  rownames(summary) <- c("mean", "sd")

  res_tot <- list(summary, res, importance)
  names(res_tot) <- c("summary","confusion", "importance")
  return(res_tot)
}

