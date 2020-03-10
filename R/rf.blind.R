#' Blind cross validation for random forest
#'
#' Grows multiple random forests with blind cross validation: the algorithm is trained
#' a specific part of the dataset, and prediction are done on another part of the dataset.
#'
#' @param tab An abundance table containing samples in columns and OTUs/ASV in rows.
#' @param treat A vector containing the class identity of each sample.
#' @param train.id A charecter sting to be searched in samples names that will be used for training.
#' @param mtry The mtry parameter to be passed to the \code{ranger} function.
#' See \code{ranger} documentation for details.
#' @param n_tree The number of tree to grow. The default is \code{500}.
#' @param n_forest The number of forests to grow. The default is \code{100}.
#' @param seed A number to set seed before sampling samples in the n-folding process
#' and before growing each forest. The default is \code{NA}.
#'
#'
#' @return Returns a list object containing the confusion matrix, the error rate, the sensitivity
#' the precision as well as the variable importance obtained for each of the \code{n_forest} grown
#' forests.
#'
#' @import ranger

# 2020-02-27
# Marine C. Cambon

rf.blind <- function(tab, treat,
                     train.id = NA,
                     mtry = NULL,
                     n_tree = 500,
                     n_forest = 100) {
  train.idx <- grep(train.id, colnames(tab))
  tab <- data.frame("treat" = treat, t(tab))
  train <- tab[train.idx, ]
  test <- tab[-train.idx, ]
  pred.irri <- error <- rate <- NULL
  res <- data.frame()
  importance <- list()
  for (i in 1:n_forest) {
    message("Growing forest number ", i, "...")
    #set.seed(140)
    rg.irri <- ranger::ranger(treat ~ ., data = train,
                      num.trees = n_tree,
                      mtry = mtry,
                      importance = "impurity")

    pred.irri <- stats::predict(rg.irri, data = test)
    error <- data.frame(table(pred.irri$predictions, test$treat))
    err_rate <- sum(test$treat != pred.irri$predictions)/nrow(test)
    TN <- error[error$Var1=="irr" & error$Var2=="irr","Freq"]
    TP <- error[error$Var1=="non-irr" & error$Var2=="non-irr","Freq"]
    FN <- error[error$Var1=="non-irr" & error$Var2=="irr","Freq"]
    FP <- error[error$Var1=="irr" & error$Var2=="non-irr","Freq"]
    sensitivity <- TP/(TP+FN)
    precision <- TP/(TP+FP)
    res <- rbind(res, cbind(TP, TN, FP, FN, err_rate, sensitivity, precision))
    importance[[i]] <- rg.irri$variable.importance
  }
  res_tot <- list(res, importance)
  names(res_tot) <- c("confusion", "importance")
  return(res_tot)
}
