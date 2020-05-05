#' k-fold cross validation for random forest
#'
#' Splits the dataset in \code{k} and grows \code{k} random forests for classification, using alternatively each
#' of the \code{k} parts of the dataset to make predictions, while the other \code{k-1} parts are used for the training.
#'
#' @param tab An abundance table containing samples in columns and OTUs/ASV in rows.
#' @param treat A boolean vector containing the class identity of each sample, i.e. the treatment to predict.
#' This means that you should pick a class as a reference for the calculation of precision and sensitivity.
#' @param k.fold A number of fold to be applied for k-fold cross-valisation.
#' @param mtry The mtry parameter to be passed to the \code{ranger} function. See \code{ranger} documentation for details.
#' @param n.tree The number of tree to grow in each of the \code{k} forests. The default is \code{500}.
#' @param seed A number to set the seed before before growing each forest. The default is \code{NULL}.
#'
#' @return A list object containing:
#' \itemize{
#'   \item a summary table with the number of true positives (TP), true negatives (TN), false positives (FP) and false negatives (FN)
#' the error rate, the sensistivity \eqn{TP/(TP + FN)}, and the precision \eqn{TP/(TP + FP)}
#'   \item The confusion matrix
#'   \item \code{n.forest} tables containing Gini index for each variable in each of the \code{n.forest} grown forests.
#'   This index gives the variable importance for classification.
#' }
#'
#' @examples
#' # Coming soon!
#'
#' @import ranger
#' @export rf.kfold

# 2020-02-27
# Marine C. Cambon


rf.kfold <- function(tab, treat,
                     k.fold = 5,
                     mtry = NULL,
                     n.tree = 500,
                     seed = NULL) {
  if(any(!treat %in% c("TRUE", "FALSE"))) stop("treat is not a boolean vector")
  treat <- ifelse(treat, "positive", "negative")
  # Preparing training IDs and dataframe
  if (!is.na(seed)) set.seed(seed)
  train.idx <- sample(rep(1:k.fold, 1/k.fold * ncol(tab)), replace = F)
  tab_agg <- data.frame("treat" = treat, t(tab))

  res <- NULL
  importance <- list()
  err_mean <- err_sd <- NULL
  TP_mean <- TN_mean <- FP_mean <- FN_mean <- NULL
  TP_sd <- TN_sd <- FP_sd <- FN_sd <- NULL
  sensitivity_mean <- precision_mean <- NULL
  sensitivity_sd <- precision_sd <- NULL

  for (i in 1:k.fold) {
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
    TN <- tmp[tmp$Var1=="negative" & tmp$Var2=="negative","Freq"]
    TP <- tmp[tmp$Var1=="positive" & tmp$Var2=="positive","Freq"]
    FN <- tmp[tmp$Var1=="positive" & tmp$Var2=="negative","Freq"]
    FP <- tmp[tmp$Var1=="negative" & tmp$Var2=="positive","Freq"]
    error <- sum(test$treat != pred$predictions)/nrow(test)
    sensitivity <- TP/(TP+FN)
    precision <- TP/(TP+FP)
    res <- rbind(res,c(TN,TP,FN,FP,error,sensitivity,precision))
    importance[[i]] <- ranger::importance(rg)
  }
  colnames(res) <- c("TN","TP","FN","FP","error","sensitivity","precision")
  rownames(res) <- paste("kfold_", 1:k.fold,sep="")
  names(importance) <- paste("kfold_", 1:k.fold,sep="")

  summary <- rbind(apply(res,2,mean),apply(res,2,sd))
  rownames(summary) <- c("mean", "sd")

  res_tot <- list(summary, res, importance)
  names(res_tot) <- c("summary","confusion", "importance")
  return(res_tot)
}

