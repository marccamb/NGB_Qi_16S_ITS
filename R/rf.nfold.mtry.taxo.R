# Funtion rf.mtry.taxo() 
# runs random forest classification with several 
#  taxonomic level and mtry parameters and returns
#  error rate obtained by n-foldinf cross-validation
#
# 2020-02-27
# Marine C. Cambon

rf.nfold.mtry.taxo <- function(tab, tax.table, treat,
                         n_fold = 5, 
                         n_mtry = 5,
                         train.id = NA, 
                         seed=1409) {
  # # For debugging:
  # l <- "ASV"
  # i <- 1
  # n <- 1
  # n_fold <- 5
  # n_mtry <- 5
  # seed <- 1409
  # tax.table <- taxo
  # variable <- "irrigation"
  # treat <- z$irrigation
  # train.id <- "-M-"
  
  # define train samples indices
  if (is.na(train.id)) {
    set.seed(seed)
    train.idx <- sample(rep(1:n_fold, 1/n_fold * ncol(tab)), replace = F)
  } else {
    train.idx <- grep(train.id, colnames(tab))
  }
  
  res_tot <- list()
  for (l in c("ASV", "genus", "family", "order", "class")) {
    if(l=="ASV") {
      tab_agg <- tab
    } else {
      tab_agg <- agg.table.taxo(tab, tax.lvl = l, taxo)
    }
    tab_agg <- data.frame("treat" = treat, t(tab_agg))
    if (n_mtry+1>ncol(tab_agg)) n_mtry <- ncol(tab_agg)-1
    
    # Defining all the variables that will be returned for each
    #   taxo level
    mtry <- 1:n_mtry*(ncol(tab_agg)-1)/n_mtry
    for (n in 1:n_mtry) {
      # Defining the variables for a given mtry
      err_mean <- NULL
      err_sd <- NULL
      TP_mean <- TN_mean <- FP_mean <- FN_mean <- NULL
      TP_sd <- TN_sd <- FP_sd <- FN_sd <- NULL
      sensitivity_mean <- precision_mean <- NULL
      sensitivity_sd <- precision_sd <- NULL
      confusion <- data.frame()
      rate <- sensitivity <- precision <- NULL
      # if (is.na(train.id)) {
        for (i in 1:n_fold) {
          train <- tab_agg[train.idx != i, ]
          test <- tab_agg[train.idx == i, ]
          set.seed(seed)
          rg.irri <- ranger(treat ~ ., data = train, 
                            num.trees = 500,
                            mtry = function(x) n*x/n_mtry,
                            importance = "impurity")
          
          pred.irri <- predict(rg.irri, data = test)
          error <- data.frame(table(pred.irri$predictions, test$treat))
          TN <- error[error$Var1=="irr" & error$Var2=="irr","Freq"]
          TP <- error[error$Var1=="non-irr" & error$Var2=="non-irr","Freq"]
          FN <- error[error$Var1=="non-irr" & error$Var2=="irr","Freq"]
          FP <- error[error$Var1=="irr" & error$Var2=="non-irr","Freq"]
          confusion <- rbind(confusion, error)
          rate <- c(rate, (FP+FN)/(FP+FN+TP+TN))
          sensitivity <- c(sensitivity, TP/(TP+FN))
          precision <- c(precision, TP/(TP+FP))
        }
      # } else {
      #   train <- tab_agg[train.idx, ]
      #   test <- tab_agg[-train.idx, ]
      #   set.seed(seed)
      #   rg.irri <- ranger(treat ~ ., data = train, 
      #                     num.trees = 500,
      #                     mtry = function(x) n*x/n_mtry,
      #                     importance = "impurity")
      #   
      #   pred.irri <- predict(rg.irri, data = test)
      #   error <- data.frame(table(pred.irri$predictions, test$treat))
      #   TN <- error[error$Var1=="irr" & error$Var2=="irr","Freq"]
      #   TP <- error[error$Var1=="non-irr" & error$Var2=="non-irr","Freq"]
      #   FN <- error[error$Var1=="non-irr" & error$Var2=="irr","Freq"]
      #   FP <- error[error$Var1=="irr" & error$Var2=="non-irr","Freq"]
      #   confusion <- rbind(confusion, error)
      #   rate <- c(rate, (FP+FN)/(FP+FN+VP+VN))
      #   sensitivity <- c(sensitivity, TP/(TP+FN))
      #   precision <- c(precision, TP/(TP+FP))
      # }
      TN_mean <- c(TN_mean, mean(confusion[confusion$Var1=="irr" & confusion$Var2=="irr","Freq"]))
      TP_mean <- c(TP_mean, mean(confusion[confusion$Var1=="non-irr" & confusion$Var2=="non-irr","Freq"]))
      FN_mean <- c(FN_mean, mean(confusion[confusion$Var1=="non-irr" & confusion$Var2=="irr","Freq"]))
      FP_mean <- c(FP_mean, mean(confusion[confusion$Var1=="irr" & confusion$Var2=="non-irr","Freq"]))
      TN_sd <- c(TN_sd, sd(confusion[confusion$Var1=="irr" & confusion$Var2=="irr","Freq"]))
      TP_sd <- c(TP_sd, sd(confusion[confusion$Var1=="non-irr" & confusion$Var2=="non-irr","Freq"]))
      FN_sd <- c(FN_sd, sd(confusion[confusion$Var1=="non-irr" & confusion$Var2=="irr","Freq"]))
      FP_sd <- c(FP_sd, sd(confusion[confusion$Var1=="irr" & confusion$Var2=="non-irr","Freq"]))
      err_mean <- c(err_mean, mean(rate))
      err_sd <- c(err_sd, sd(rate))
      sensitivity_mean <- c(sensitivity_mean, mean(sensitivity))
      precision_mean <- c(precision_mean, mean(precision))
      sensitivity_sd <- c(sensitivity_sd, sd(sensitivity))
      precision_sd <- c(precision_sd, sd(precision))
    }
    res_tot[[l]] <- data.frame(cbind(TP_mean, TP_sd, 
                                     TN_mean, TN_sd, 
                                     FP_mean, FP_sd, 
                                     FN_mean, FN_sd, 
                                     err_mean, err_sd,
                                     sensitivity_mean,
                                     sensitivity_sd,
                                     precision_mean,
                                     precision_sd,
                                     mtry))
    message(l, " lvl is done\n")
  }
  return(res_tot)
}

rf.nfold <- function(tab, treat,
                     n_fold = 5, 
                     mtry = NULL,
                     seed=1409) {
  set.seed(seed)
  train.idx <- sample(rep(1:n_fold, 1/n_fold * ncol(tab)), replace = F)
  tab <- data.frame("treat" = treat, t(tab))  
  pred.irri <- error <- rate <- NULL
  res <- data.frame()
  importance <- list()
  for (i in 1:n_fold) {
    train <- tab[train.idx != i, ]
    test <- tab[train.idx == i, ]
    set.seed(seed)
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

rf.define.training <- function(tab, treat,
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