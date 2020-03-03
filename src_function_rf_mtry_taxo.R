# Funtion rf.mtry.taxo() 
# runs random forest classification with several 
#  taxonomic level and mtry parameters and returns
#  error rate obtained by n-foldinf cross-validation
#
# 2020-02-27
# Marine C. Cambon

rf.mtry.taxo <- function(tab, tax.table, treat,
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
    mtry <- 1:n_mtry*(ncol(tab_agg)-1)/n_mtry
    err_mean <- NULL
    err_sd <- NULL
    TP <- TN <- FP <- FN <- NULL
    
    for (n in 1:n_mtry) {
      error <- data.frame()
      rate <- NULL
      if (is.na(train.id)) {
        for (i in 1:n_fold) {
          train <- tab_agg[train.idx != i, ]
          test <- tab_agg[train.idx == i, ]
          set.seed(1409)
          rg.irri <- ranger(treat ~ ., data = train, 
                            num.trees = 500,
                            mtry = function(x) n*x/n_mtry,
                            importance = "impurity")
          
          pred.irri <- predict(rg.irri, data = test)
          #data.frame(table(pred.irri$predictions, test$treat))
          error <- rbind(error, data.frame(table(pred.irri$predictions, test$treat)))
          rate <- c(rate, sum(test$treat != pred.irri$predictions)/nrow(test))
          #print(table(test$irrigation, pred.irri$predictions))
        }
      } else {
        train <- tab_agg[train.idx, ]
        test <- tab_agg[-train.idx, ]
        set.seed(1409)
        rg.irri <- ranger(treat ~ ., data = train, 
                          num.trees = 500,
                          mtry = function(x) n*x/n_mtry,
                          importance = "impurity")
        
        pred.irri <- predict(rg.irri, data = test)
        data.frame(table(pred.irri$predictions, test$treat))
        error <- rbind(error, data.frame(table(pred.irri$predictions, test$treat)))
        rate <- c(rate, sum(test$treat != pred.irri$predictions)/nrow(test))
      }
      TN <- c(TN, mean(error[error$Var1=="irr" & error$Var2=="irr","Freq"]))
      TP <- c(TP, mean(error[error$Var1=="non-irr" & error$Var2=="non-irr","Freq"]))
      FN <- c(FN, mean(error[error$Var1=="non-irr" & error$Var2=="irr","Freq"]))
      FP <- c(FP, mean(error[error$Var1=="irr" & error$Var2=="non-irr","Freq"]))
      err_mean <- c(err_mean, mean(rate))
      err_sd <- c(err_sd, sd(rate))
    }
    res_tot[[l]] <- data.frame(cbind(TP, TN, FP, FN, err_mean, err_sd, mtry))
    res_tot[[l]]$sensitivity <- with(res_tot[[l]], TP/(TP+FN))
    res_tot[[l]]$precision <- with(res_tot[[l]], TP/(TP+FP))
    message(l, " lvl is done\n")
  }
  return(res_tot)
}