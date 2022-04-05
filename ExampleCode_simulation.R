load("OneSimulationData.RData")
source("functions_existings.R")
library(CAMLU)
library(keras)
library(mclust)

### run CAMLU
label_01 <- CAMLU(x_train = as.matrix(simdata$fdata_train),
                  x_test = as.matrix(simdata$fdata_test),
                  ngene=3000,lognormalize=TRUE)
acc_A <- sum(label_01 == simdata$labels_2)/length(label_01)
ARI_A <- adjustedRandIndex(label_01, simdata$labels_2)

label_all <- CAMLU(x_train = as.matrix(simdata$fdata_train),
                   x_test = as.matrix(simdata$fdata_test),
                   y_train = simdata$labels_train,
                   full_annotation = TRUE,
                   ngene=3000,
                   lognormalize=TRUE)
fulllabel <- label_all$label_full
truelabel <- simdata$labels_test
truelabel[simdata$labels_test == "cd14_monocytes"] <- "Unknown"
acc_A2 <- sum(fulllabel == truelabel)/length(truelabel)
ARI_A2 <- adjustedRandIndex(fulllabel, truelabel)


### scmap-cluster
message("scmap-cluster")
scmapcluster_out <- GetLabels(train_count = simdata$fdata_train, 
                              train_label = simdata$labels_train,
                              test_count = simdata$fdata_test,
                              method = "scmapcluster")
newscmapcluster <- scmapcluster_out
newscmapcluster[scmapcluster_out == "unassigned"] <- 1
newscmapcluster[scmapcluster_out != "unassigned"] <- 0
acc_SC <- sum(newscmapcluster == simdata$labels_2)/length(newscmapcluster)
ARI_SC <- adjustedRandIndex(newscmapcluster, simdata$labels_2)

newscmapcluster2 <- scmapcluster_out
newscmapcluster2[newscmapcluster2 == "unassigned"] = "Unknown"
acc_SC2 <- sum(newscmapcluster2 == truelabel)/length(truelabel)
ARI_SC2 <- adjustedRandIndex(newscmapcluster2, truelabel)



### scmap-cell
message("scmap-cell")
scmapcell_out <- GetLabels(train_count = simdata$fdata_train, 
                           train_label = simdata$labels_train,
                           test_count = simdata$fdata_test,
                           method = "scmapcell")
newscmapcell <- scmapcell_out
newscmapcell[scmapcell_out == "unassigned"] <- 1
newscmapcell[scmapcell_out != "unassigned"] <- 0
acc_SE <- sum(newscmapcell == simdata$labels_2)/length(newscmapcell)
ARI_SE <- adjustedRandIndex(newscmapcell, simdata$labels_2)

newscmapcell2 <- scmapcell_out
newscmapcell2[newscmapcell2 == "unassigned"] = "Unknown"
acc_SE2 <- sum(newscmapcell2 == truelabel)/length(truelabel)
ARI_SE2 <- adjustedRandIndex(newscmapcell2, truelabel)

### CHETAH
message("CHETAH")
chetah_out <- GetLabels(train_count = simdata$fdata_train, 
                        train_label = simdata$labels_train,
                        test_count = simdata$fdata_test,
                        method = "CHETAH")
newchetah <- chetah_out
newchetah[chetah_out == "unassigned"] <- 1
newchetah[chetah_out != "unassigned"] <- 0
acc_CT <- sum(newchetah == simdata$labels_2)/length(newchetah)
ARI_CT <- adjustedRandIndex(newchetah, simdata$labels_2)

newchetah2 <- chetah_out
newchetah2[newchetah2 == "unassigned"] = "Unknown"
acc_CT2 <- sum(newchetah2 == truelabel)/length(truelabel)
ARI_CT2 <- adjustedRandIndex(newchetah2, truelabel)


### CHETAH
message("scPred")
scPred_out <- GetLabels(train_count = simdata$fdata_train, 
                        train_label = simdata$labels_train,
                        test_count = simdata$fdata_test,
                        method = "scPred")
newscPred <- scPred_out
newscPred[scPred_out == "unassigned"] <- 1
newscPred[scPred_out != "unassigned"] <- 0
acc_SP <- sum(newscPred == simdata$labels_2)/length(newscPred)
ARI_SP <- adjustedRandIndex(newscPred, simdata$labels_2)

tttable <- as.matrix(table(newscPred, simdata$labels_2))
Precision_SP <- diag(tttable)/rowSums(tttable)
Recall_SP <- diag(tttable)/colSums(tttable)

newscPred2 <- scPred_out
newscPred2[newscPred2 == "unassigned"] = "Unknown"
acc_SP2 <- sum(newscPred2 == truelabel)/length(truelabel)
ARI_SP2 <- adjustedRandIndex(newscPred2, truelabel)


### scPred
message("scPred")
scPred_out <- GetLabels(train_count = simdata$fdata_train, 
                        train_label = simdata$labels_train,
                        test_count = simdata$fdata_test,
                        method = "scPred")
newscPred <- scPred_out
newscPred[scPred_out == "unassigned"] <- 1
newscPred[scPred_out != "unassigned"] <- 0
acc_SP <- sum(newscPred == simdata$labels_2)/length(newscPred)
ARI_SP <- adjustedRandIndex(newscPred, simdata$labels_2)

newscPred2 <- scPred_out
newscPred2[newscPred2 == "unassigned"] = "Unknown"
acc_SP2 <- sum(newscPred2 == truelabel)/length(truelabel)
ARI_SP2 <- adjustedRandIndex(newscPred2, truelabel)


allres <- matrix(NA, 5, 4)
allres[1, ] <- c(acc_A, ARI_A, acc_A2, ARI_A2)
allres[2, ] <- c(acc_SC, ARI_SC, acc_SC2, ARI_SC2)
allres[3, ] <- c(acc_SE, ARI_SE, acc_SE2, ARI_SE2)
allres[4, ] <- c(acc_CT, ARI_CT, acc_CT2, ARI_CT2)
allres[5, ] <- c(acc_SP, ARI_SP, acc_SP2, ARI_SP2)

rownames(allres) <- c("our", "scmap-cluster",
                      "scmap-cell", "CHETAH", "scPred")
colnames(allres) <- c("KnownUnkonwn_Accuracy",
                      "KnownUnknown_ARI",
                      "Fulllabel_Accuracy",
                      "Fulllabel_ARI")
allres

