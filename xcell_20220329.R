libarary(xCell)
library(Seurat)
expdata=expr[['Spatial']]@data
data=xCellAnalysis(expdata,rnaseq=F,parallel.sz = 16,save.raw=TRUE,file.name='test')
expr1=expr
for(i in rownames(data)[1:64]){
  expr1[[i]]=data[i,colnames(expr1)]
}
colnames(expr1@meta.data)[10:73]=rownames(data)[1:64]
#for(i in rownames(data))
###"X3.18S08906.002" "X4.19S58966.002"
slices=names(expr1@images)
for(slice in slices){
  p1=SpatialFeaturePlot(expr1,features=colnames(expr1@meta.data)[10:73],images=slice,
                      ncol=5,alpha = c(0.1,1),crop=FALSE)
  ggsave(paste('xcell_',slice,'.png',sep=""),p1,width=16,height=70,limitsize = FALSE)
  ggsave(paste('xcell_',slice,'.pdf',sep=""),p1,width=16,height=70,limitsize = FALSE)
}
expr2=expr
for(i in rownames(rawgsvadata)[1:64]){
  expr2[[i]]=as.numeric(rawgsvadata[i,colnames(expr2)][1,])
}
colnames(expr2@meta.data)[10:73]=rownames(rawgsvadata)[1:64]
#for(i in rownames(data))
###"X3.18S08906.002" "X4.19S58966.002"
slices=names(expr2@images)
for(slice in slices){
  p1=SpatialFeaturePlot(expr2,features=colnames(expr2@meta.data)[10:73],images=slice,
                        ncol=5,alpha = c(0.1,1),crop=FALSE)
  ggsave(paste('gava_',slice,'.png',sep=""),p1,width=16,height=70,limitsize = FALSE)
  ggsave(paste('gsva_',slice,'.pdf',sep=""),p1,width=16,height=70,limitsize = FALSE)
}
geneSets=signatures
rawgsvadata=rawgsvaAnalysis(expdata,geneSets,parallel.sz=20,parallel.type='SOCK',file.name='rawGSVA')

####################20220330 msigdb  xcell
msigdbxcellscores=xCellAnalysis(expdata,signatures = msigdb_used, 
                   rnaseq=F,parallel.sz = 16,save.raw=TRUE,file.name='msigdbxcellscores')
spatialplot(expr,filename='msigdb_xcell',scoredata=msigdbxcellscores)
  
####################20220330 msigdb gsva
# scores <- rawEnrichmentAnalysis(expr, signatures, genes, fn, parallel.sz = parallel.sz, parallel.type = "SOCK")
msigdbgsvadata=rawgsvaAnalysis(expdata,msigdb_used,parallel.sz=20,parallel.type='SOCK',file.name='msigdbGSVA')
spatialplot(expr,filename='msigdb_gsva',scoredata=msigdbgsvadata)

spatialplot<-function(expr,filename,scoredata){
  expr1=expr
  print(head(colnames(expr1)))
  for(i in rownames(scoredata)[1:nrow(scoredata)]){
    expr1[[i]]=scoredata[i,colnames(expr1)]
  }
  colnames(expr1@meta.data)[(ncol(expr@meta.data)+1):ncol(expr1@meta.data)]=rownames(scoredata)[1:nrow(scoredata)]
  #for(i in rownames(data))
  ###"X3.18S08906.002" "X4.19S58966.002"
  slices=names(expr1@images)
  for(slice in slices){
    print(slice)
    #[10:ncol(expr1@meta.data)]
    p1=SpatialPlot(expr1,features=colnames(expr1@meta.data)[10:ncol(expr1@meta.data)],images=slice,
                          ncol=3,alpha = c(0.1,1),crop=FALSE,pt.size.factor = 1.3,label.size=10)
    height=5.6*(nrow(scoredata)/3)
    ggsave(paste(filename,'_',slice,'.png',sep=""),p1,width=16,height=height,limitsize = FALSE)
    ggsave(paste(filename,'_',slice,'.pdf',sep=""),p1,width=16,height=height,limitsize = FALSE)
  }
}
rawgsvaAnalysis<-function(expdata,geneSets,parallel.sz=30,parallel.type='SOCK',file.name='rawGSVA'){
  expr <- apply(expdata, 2, rank)
  #gs <- gsva(expr, method = "ssgsea", geneSets, mx.diff = TRUE,parallel.sz=20)
  scores <- GSVA::gsva(expr, geneSets, method = "ssgsea",mx.diff = TRUE,parallel.sz = parallel.sz)
  scores = scores - apply(scores, 1, min)  ##加上这一行命令脑子是秀逗了吗
  if('%'%in%rownames(scores)[1]){
    cell_types <- unlist(strsplit(rownames(scores), "%"))
    cell_types <- cell_types[seq(1, length(cell_types), 3)]
    agg <- aggregate(scores ~ cell_types, FUN = mean)
    rownames(agg) <- agg[, 1]
    scores <- agg[, -1]
  }
  write.table(scores, file = file.name, sep = "\t", col.names = NA,
              quote = FALSE)
  return(scores)
}


xCellAnalysis<-function (expr, signatures = NULL, genes = NULL, spill = NULL, 
          rnaseq = TRUE, file.name = NULL, scale = TRUE, alpha = 0.5, 
          save.raw = FALSE, parallel.sz = 4, parallel.type = "SOCK", 
          cell.types.use = NULL){
  if (is.null(signatures)) 
    signatures = xCell.data$signatures
  if (is.null(genes)) 
    genes = xCell.data$genes
  if (is.null(spill)) {
    if (rnaseq == TRUE) {
      spill = xCell.data$spill
    }
    else {
      spill = xCell.data$spill.array
    }
  }
  if (is.null(file.name) || save.raw == FALSE) {
    fn <- NULL
  }
  else {
    fn <- paste0(file.name, "_RAW.txt")
  }
  if (!is.null(cell.types.use)) {
    A = intersect(cell.types.use, rownames(spill$K))
    if (length(A) < length(cell.types.use)) {
      return("ERROR - not all cell types listed are available")
    }
  }
  scores <- rawEnrichmentAnalysis(expr, signatures, genes, 
                                  fn, parallel.sz = parallel.sz, parallel.type = "SOCK")
  scores.transformed <- transformScores(scores, spill$fv, scale)
  if (is.null(file.name)) {
    fn <- NULL
  }
  else {
    fn <- file.name
  }
  if (is.null(cell.types.use)) {
    scores.adjusted <- spillOver(scores.transformed, spill$K, alpha, fn)
    scores.adjusted = microenvironmentScores(scores.adjusted)
  }
  else {
    scores.adjusted <- spillOver(scores.transformed[cell.types.use, 
    ], spill$K, alpha, fn)
  }
  return(scores.adjusted)
}
#scores.adjusted <- spillOver(scores.transformed, spill$K, alpha, fn)
spill = xCell.data$spill.array
spillOver<-function (transformedScores, K, alpha = 0.5, file.name = NULL){
  K <- K * alpha
  diag(K) <- 1
  rows <- rownames(transformedScores)[rownames(transformedScores) %in% rownames(K)]
  scores <- apply(transformedScores[rows, ], 2, function(x) pracma::lsqlincon(K[rows, rows], x, lb = 0))
  scores[scores < 0] = 0
  rownames(scores) <- rows
  if (!is.null(file.name)) {
    scores = round(scores * 10000)
    scores = scores/10000
    write.table(scores, file = file.name, sep = "\t", col.names = NA, 
                quote = FALSE)
  }
  return(scores)
}
microenvironmentScores<-function (adjustedScores){
  ImmuneScore = apply(adjustedScores[c("B-cells", "CD4+ T-cells", 
                                       "CD8+ T-cells", "DC", "Eosinophils", "Macrophages", "Monocytes", 
                                       "Mast cells", "Neutrophils", "NK cells"), ], 2, sum)/1.5
  StromaScore = apply(adjustedScores[c("Adipocytes", "Endothelial cells", 
                                       "Fibroblasts"), ], 2, sum)/2
  MicroenvironmentScore = ImmuneScore + StromaScore
  adjustedScores = rbind(adjustedScores, ImmuneScore, StromaScore, 
                         MicroenvironmentScore)
}
rawEnrichmentAnalysis<-function (expr, signatures, genes, file.name = NULL, parallel.sz = 4, 
          parallel.type = "SOCK"){
  shared.genes <- intersect(rownames(expr), genes)
  print(paste("Num. of genes:", length(shared.genes)))
  expr <- expr[shared.genes, ]
  if (dim(expr)[1] < 5000) {
    print(paste("ERROR: not enough genes"))
    return - 1
  }
  expr <- apply(expr, 2, rank)
  scores <- GSVA::gsva(expr, signatures, method = "ssgsea", 
                       ssgsea.norm = FALSE, parallel.sz = parallel.sz, parallel.type = parallel.type)
  scores = scores - apply(scores, 1, min)
  saveRDS(scores,'xcellmsigdbscores.RDS')
  if('%'%in%rownames(scores)[1]){
    cell_types <- unlist(strsplit(rownames(scores), "%"))
    cell_types <- cell_types[seq(1, length(cell_types), 3)]
    agg <- aggregate(scores ~ cell_types, FUN = mean)
    rownames(agg) <- agg[, 1]
    scores <- agg[, -1]
  }else{
    scores=scores
  }
  if (!is.null(file.name)) {
    write.table(scores, file = file.name, sep = "\t", col.names = NA, 
                quote = FALSE)
  }
  scores
}

####################GSVA分析信息

geneSets=lapply(signatures,function(t)t@geneIds)
names(geneSets)=names(signatures)
#########BHT201002 查找spatial文件夹
# agg=read.delim('/data7/guohua/Visium/BHT201002_40samples/RNA/AGG_libraries.csv',sep=",",as.is=T,
#               header=T,check.names=F)
