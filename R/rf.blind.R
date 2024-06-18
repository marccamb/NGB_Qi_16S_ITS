#' Non-random cross validation for random forest
#'
#' Grows multiple random forests with non-random cross validation: the algorithm is trained on
#' a specific part of the dataset, and predictions are done on another part of the dataset.
#'
#' @param tab An abundance or presence absence table containing samples in columns and OTUs/ASV in rows.
#' @param treat A boolean vector containing the class identity of each sample, i.e. the treatment to predict.
#' This means that you should pick a class as a reference for the calculation of precision and sensitivity.
#' @param test.id A charecter sting to be searched in samples names that will be used for testing.
#' Can be a regular expression. Can alernatively be a boolean vector saying wether or not each sample
#' is part of the testing or training dataset (TRUE for testing samples, FALSE for training samples), or a character
#' vector containing the testing sample names.
#' @param mtry The mtry parameter to be passed to the \code{ranger} function.
#' See \code{ranger} documentation for details.
#' @param n.tree The number of tree to grow. The default is \code{500}.
#' @param n.forest The number of forests to grow. The default is \code{10}.
#' @param importance_p A boolean defining if the p-value should be computed for the importance of
#' variable. For now, the importance is the Gini index, and the p-value is estimated by permutation with
#' the Altmann method. See ranger documentation for details
#' @param seed A number to set the seed before growing the forest. Only meaningful
#' if n.forest == 1. The default is \code{NULL}.
#'
#'@return A list object containing:
#' \itemize{
#'   \item a summary table with the number of true positives (TP), true negatives (TN), false positives (FP) and false negatives (FN)
#' the error rate, the sensistivity \eqn{TP/(TP + FN)}, and the precision \eqn{TP/(TP + FP)}
#'   \item The confusion matrix
#'   \item \code{n.forest} tables containing Gini index for each variable in each of the \code{n.forest} grown forests.
#'   This index gives the variable importance for classification.
#' }
#'
#' @import ranger
#' @export rf.blind

# 2020-02-27
# Marine C. Cambon

rf.blind <- function(tab, treat,
                     test.id,
                     mtry = NULL,
                     n.tree = 500,
                     n.forest = 10,
                     importance_p = F,
                     seed=NULL) {
  if(class(treat) != "logical") stop("treat is not a boolean vector")
  treat <- ifelse(treat, "positive", "negative")
  treat <- as.factor(treat)

  if(length(test.id)==1) {
    test.idx <- grep(test.id, colnames(tab))
  } else {
    if(class(test.id) == "logical") {
      test.idx <- which(test.id)
    } else {
      test.idx <- which(colnames(tab) %in% test.id)
    }
  }
  if(length(test.idx)==1) warning("The testing dataset only contains 1 sample")
  if(length(test.idx)==0) stop("test.id does not match sample names")

  tab <- data.frame("treat" = treat, t(tab))
  train <- tab[-test.idx, ]
  test <- tab[test.idx, ]
  res <- data.frame()
  importance <- list()
  message("Growing ", n.forest, " forests...")
  for (i in 1:n.forest) {
    if(n.forest == 1) set.seed(seed)
    rg <- ranger::ranger(treat ~ ., data = train,
                      num.trees = n.tree,
                      mtry = mtry,
                      importance = "impurity")

    pred <- stats::predict(rg, data = test)
    tmp <- data.frame(table(pred$predictions, test$treat))
    TN <- tmp[tmp$Var1=="negative" & tmp$Var2=="negative","Freq"]
    TP <- tmp[tmp$Var1=="positive" & tmp$Var2=="positive","Freq"]
    FN <- tmp[tmp$Var1=="positive" & tmp$Var2=="negative","Freq"]
    FP <- tmp[tmp$Var1=="negative" & tmp$Var2=="positive","Freq"]
    error <- sum(test$treat != pred$predictions)/nrow(test)
    sensitivity <- TP/(TP+FN)
    precision <- TP/(TP+FP)
    res <- rbind(res, c(rg$mtry, TP, TN, FP, FN, error, sensitivity, precision))
    if (importance_p) {importance[[i]] <- ranger::importance_pvalues(rg, method = "altmann",
                                                  formula=treat ~ .,
                                                  data = train)
    } else {
      importance[[i]] <- rg$variable.importance
    }
  }
  colnames(res) <- c("mtry","TP","TN","FP","FN","error","sensitivity","precision")
  message("Done!")
  summary <- rbind(apply(res,2,mean),apply(res,2,sd))
  rownames(summary) <- c("mean", "sd")

  res_tot <- list(summary, res, importance)
  names(res_tot) <- c("summary", "confusion", "importance")
  return(res_tot)
}
