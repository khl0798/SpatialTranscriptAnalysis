#######20220811样本配对分析
# 017S01321-005和520S22781-003配对，617S33388-001和319S04687-002配对。分析的话就还是按照之前的流程即可。
# S13-PT与S13-MT配对，S16-PT与S16-MT配对
library(Seurat)
runScHarmony<-function(sampleName, dataPath, savePath, hvgGenesCombine=FALSE){
  library("hdf5r")
  imgpath <- paste(datapath, sampleName, "/outs/spatial", sep = "")
  exprf <- paste(datapath, sampleName, "/outs/filtered_feature_bc_matrix.h5", 
                 sep = "")
  expr.data <- Seurat::Read10X_h5(filename = exprf)
  expr <- Seurat::CreateSeuratObject(counts = expr.data, project = sampleName, 
                                     assay = "Spatial")
  expr$slice <- 1
  expr$region <- sampleName
  img <- Seurat::Read10X_Image(image.dir = imgpath)
  Seurat::DefaultAssay(object = img) <- "Spatial"
  img <- img[colnames(x = expr)]
  expr[[sampleName]] <- img
  plot1 <- VlnPlot(expr, features = "nCount_Spatial", pt.size = 0.1) + 
    NoLegend()
  plot2 <- SpatialFeaturePlot(expr, features = "nCount_Spatial") + 
    theme(legend.position = "right")
  pdf(paste(sampleName, "_Total_UMI_in_spots.pdf", sep = ""), 
      width = 8, height = 4, onefile = F)
  print(plot_grid(plot1, plot2))
  dev.off()
  plot1 <- VlnPlot(expr, features = "nFeature_Spatial", pt.size = 0.1) + 
    NoLegend()
  plot2 <- SpatialFeaturePlot(expr, features = "nFeature_Spatial") + 
    theme(legend.position = "right")
  pdf(paste(sampleName, "_Total_genes_in_spots.pdf", sep = ""), 
      width = 8, height = 4, onefile = F)
  print(plot_grid(plot1, plot2))
  dev.off()
  # expr <- SCTransform(expr, assay = "Spatial", return.only.var.genes = FALSE, 
                      # verbose = FALSE)
  expr <- NormalizeData(expr,normalization.method = 'LogNormalize',scale.factor=10000)
  expr <- FindVariableFeatures(expr,selection.method = "vst", nfeatures = 2000, verbose = F)
  expr <- ScaleData(expr,vars.to.regress = NULL,verbose = F)
  saveRDS(expr, file = file.path(savePath, paste(sampleName, ".expr.RDS", sep = "")))
}
##########runharmony method
runSpCombination1<-function (sampleNames, savePath, pc.use = 30){
  message("[", Sys.time(), "] -----: sample data combination")
  expr.list <- list()
  sample.ident <- c()
  sample.Var <- c()
  for (i in 1:length(sampleNames)) {
    sampleName <- sampleNames[i]
    print(sampleName)
    expr.list[[sampleName]] <- readRDS(paste0(savePath, "/", 
                                              sampleName, ".expr.RDS"))
    sample.ident <- c(sample.ident, rep(sampleName, dim(expr.list[[sampleName]])[2]))
    sample.Var <- c(sample.Var, VariableFeatures(expr.list[[sampleName]]))
  }
  sample.ident <- as.factor(sample.ident)
  message("[", Sys.time(), "] -----: combine raw matrix data")
  suppressWarnings(expr <- merge(expr.list[[1]], expr.list[2:length(expr.list)]))
  # DefaultAssay(expr) <- "SCT"
  expr[["sample.ident"]] <- sample.ident
  # VariableFeatures(expr) <- unique(sample.Var)
  expr <- FindVariableFeatures(expr,selection.method = "vst", nfeatures = 2000, verbose = F)
  expr <- ScaleData(expr,vars.to.regress = NULL,verbose = F)
  expr <- RunPCA(expr, verbose = FALSE)
  message("[", Sys.time(), "] -----: combine data by harmony")
  library(harmony)
  DefaultAssay(expr)='Spatial'
  expr <- RunHarmony(expr, group.by.vars = "orig.ident",assay.use = "Spatial")
  expr <- FindNeighbors(expr, dims = 1:pc.use,reduction='harmony')
  expr <- RunUMAP(expr, dims = 1:pc.use,reduction='harmony')
  expr <- RunTSNE(expr, dims = 1:pc.use,reduction='harmony')
  expr <- run_cluster(expr, sampleNames, savePath, resolution = 0.5, 
                      clusterStashName = "default")
  return(expr)
}
