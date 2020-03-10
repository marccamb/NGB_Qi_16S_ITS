rf.blind <- function(tab, treat,
                               train.id = NA, 
                               mtry = NULL,
                               n_forest=100) {
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
    rg.irri <- ranger(treat ~ ., data = train, 
                      num.trees = 500,
                      mtry = mtry,
                      importance = "impurity")
    
    pred.irri <- predict(rg.irri, data = test)
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
