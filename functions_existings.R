### Existing method

GetLabels = function(train_count = NULL, 
                     train_label = NULL,
                     test_count,
                     nfeature_svm = 1000,
                     scPredMethod = "svmRadial",
                     method = c("scmapcluster", 
                                "scmapcell", 
                                "CHETAH", 
                                "SVM", 
                                "Seurat", 
                                "SC3",
                                "SingleR",
                                "scPred"), 
                     K = NULL) {
  
  require(SingleCellExperiment)
  
  switch(method,
         scmapcluster = label.scmap(train_count, train_label,
                                    test_count, type = "cluster"),
         scmapcell = label.scmap(train_count, train_label,
                                 test_count, type = "cell"),
         CHETAH = label.CHETAH(train_count, train_label,
                               test_count),
         SVM = label.SVM(train_count, train_label,
                               test_count, nfeature_svm),
         SC3 = cluster.SC3(test_count, K),
         Seurat = cluster.Seurat(test_count, K),
         SingleR = cluster.SingleR(train_count, train_label, test_count),
         scPred = cluster.scPred(train_count, train_label, test_count,
                                 scPredMethod)
  )
}

label.scmap <- function(train_count, train_label,
                         test_count, type = "cluster") {
  require(scmap)
  
  #train
  sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(train_count)), 
                              colData = DataFrame(cell_type1 = train_label))
  logcounts(sce) <- log2(normcounts(sce) + 1)
  # use gene names as feature symbols
  rowData(sce)$feature_symbol <- rownames(sce)
  # remove features with duplicated names
  sce <- sce[!duplicated(rownames(sce)), ]
  sce <- suppressWarnings(selectFeatures(sce, suppress_plot = TRUE))
  
  #test
  sce_test <- SingleCellExperiment(list(normcounts = as.matrix(test_count)), 
                                   colData = DataFrame(cell_type1 = rep(NA, ncol(test_count))))
  logcounts(sce_test) <- log2(normcounts(sce_test) + 1)
  rowData(sce_test)$feature_symbol <- rownames(sce_test)
  
  if(type == "cluster") {
    #scmap-cluster
    sce <- indexCluster(sce)
    scmapCluster_results <- scmapCluster(projection = sce_test,
                                         index_list = list(metadata(sce)$scmap_cluster_index))
    scmapout <- scmapCluster_results$combined_labs
  } else if(type == "cell") {
    #scmap-Cell
    sce <- suppressMessages(indexCell(sce))
    scmapCell_results <- scmapCell(projection = sce_test, 
                                   index_list = list(metadata(sce)$scmap_cell_index))
    scmapCell_clusters <- scmapCell2Cluster(scmapCell_results, 
                                            list(as.character(colData(sce)$cell_type1)))
    scmapout <- scmapCell_clusters$combined_labs
  } else {
    stop("Type need to be cluster or cell!")
  }
  
  return(scmapout)
}

label.CHETAH <- function(train_count, train_label,
                          test_count){
  
  require(CHETAH)
  
  sce <- SingleCellExperiment(assays = list(counts = train_count), 
                              colData = DataFrame(celltypes = train_label))
  #test
  sce_test <- SingleCellExperiment(assays = list(counts = test_count), 
                                   colData = DataFrame(celltypes = rep(NA, ncol(test_count))))
  
  #run classifier
  CHETAHout <- CHETAHclassifier(input = sce_test, ref_cells = sce)
  outlab <- CHETAHout$celltype_CHETAH
  outlab[which(substr(outlab, 1, 4) == "Node")] <- "unassigned"
  outlab[which(outlab == "Unassigned")] <- "unassigned"
  return(outlab)
}

label.SVM <- function(train_count, train_label,
                       test_count, nfeature = 1000) {
  require(e1071)
  
  tsd <- apply(train_count, 1, sd)
  ooidx <- order(tsd, decreasing = TRUE)[1:nfeature]
  
  train_count2 <- train_count[ooidx,]
  test_count2 <- test_count[ooidx,]
  
  uniqID <- unique(train_label)
  matlab <- match(train_label, uniqID)
  # train
  svmout <- svm(t(as.matrix(train_count2)), as.factor(matlab), kernel = "linear")
  predlab <- predict(svmout, t(as.matrix(test_count2)))
  
  outlab <- uniqID[predlab]
  return(outlab)
}

cluster.Seurat <- function(counts, K) {
  require(Seurat)
  ## PCA dimention redcution
  seuset = CreateSeuratObject( counts )
  seuset = NormalizeData(object = seuset)
  seuset = FindVariableFeatures(object = seuset)
  seuset = ScaleData(object = seuset)
  seuset = RunPCA(object = seuset)
  seuset = FindNeighbors(object = seuset)
  seuset = FindClusters(object = seuset)
  return(seuset@active.ident)
}

cluster.SC3 <- function(counts, K) {
  require(SC3)
  sce = SingleCellExperiment(
    assays = list(
      counts = as.matrix(counts),
      logcounts = log2(as.matrix(counts) + 1)
    )
  )
  rowData(sce)$feature_symbol <- rownames(sce)
  sce = sce[!duplicated(rowData(sce)$feature_symbol), ]
  sce = sc3_prepare(sce, gene_filter = FALSE)
  if( missing(K) ) { ## estimate number of clusters
    sce = sc3_estimate_k(sce)
    K = metadata(sce)$sc3$k_estimation
  }
  
  sce = sc3_calc_dists(sce)
  sce = sc3_calc_transfs(sce)
  sce = sc3_kmeans(sce, ks = K)
  sce = sc3_calc_consens(sce)
  if(ncol(counts)>5000) {
    sce = sc3_run_svm(sce, ks = K)
  }
  colTb = colData(sce)[,1]
  return(list(colTb=colTb,
              sce = sce))
}

cluster.SingleR <- function(train_count, train_label,
                            test_count) {
  
  require(scuttle)
  require(SingleR)
  
  sceTrain <- SingleCellExperiment(list(counts=train_count),
                                   colData=DataFrame(celltype=as.factor(train_label)))
  sceTest <- SingleCellExperiment(list(counts=test_count[, aout == 0]))
  sceTrain <- logNormCounts(sceTrain)
  sceTest <- logNormCounts(sceTest)
  
  pred.out <- SingleR(test = sceTest,
                      ref = sceTrain,
                      labels = sceTrain$celltype)
  return(pred.out$labels)
}

cluster.scPred <- function(train_count, train_label,
                            test_count, method = "svmRadial") {
  
  require(scPred)
  require(Seurat)
  
  reference <- CreateSeuratObject(counts = train_count)
  reference[["celltype"]] <- as.factor(train_label)
  reference <- reference %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30)
  reference <- getFeatureSpace(reference, "celltype")
  reference <- trainModel(reference, model = method)
  
  query <- CreateSeuratObject(counts = test_count)
  query <- NormalizeData(query)
  query <- scPredict(query, reference)

  return(query$scpred_prediction)
}


