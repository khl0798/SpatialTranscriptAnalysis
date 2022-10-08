options(rlib_downstream_check = FALSE)
library(Seurat)
library(Matrix)
library(dplyr)
library(plyr)
library(ggplot2)
library(spacexr)
library(future)
spatial=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/new.combin.20220422.data.RDS')
reference=readRDS('/data7/konghl/data/SpatialTranscriptome/BHT201004/analysis_1-N1-1_2-N2-1/reannotation/17clusters_20220805_singlecell/pal_scRNA.rds')
sc_anno=read.delim('sccluster_annotation.txt',sep="\t",as.is=T,header=F)
new_cluster=sc_anno[,2]
names(new_cluster)=sc_anno[,1]
reference=RenameIdents(reference,new_cluster)
# sc.reference <- readRDS("/data7/guohua/10X/BHT201002_sc/220217/28sc_GBC_final.RDS")
# Idents(sc.reference) <- "celltype"
# Spatial <- readRDS("/data7/guohua/Visium/BHT201002_40samples/analysis/combin.data.RDS")
# Spatial <- RenameCells(Spatial,new.names=gsub("1_","",colnames(Spatial)))
# tichu <- readRDS("/data7/guohua/10X/BHT201002_sc/SC_VI40_0310/myRCTD.RDS")
# bc <- colnames(tichu@spatialRNA@counts)
# bc <- gsub("1_","",bc)
# Spatial <- subset(Spatial,cell=bc)
# anno <- read.table("/data7/guohua/Visium/BHT201002_40samples/analysis/cell_sample_cluster_annotation_main.xls",header=TRUE,sep="\t")
# anno$barcodes <- gsub("1_","",anno$barcodes)

coords.file <- function(expr) {
  spatial_coord <- data.frame()
  meta <- data.frame(barcode=rownames(expr[[]]),stringsAsFactors=F)
  colnames(meta) <- "barcodeID"
  # for (i in unique(expr[["sample.ident"]])[,1]){
  for(i in names(expr@images)){
    tmp <- data.frame(expr@images[[i]]@coordinates) %>%
      tibble::rownames_to_column("barcodeID") %>%
      dplyr::mutate(imagerow_scaled = imagerow * expr@images[[i]]@scale.factors$lowres,
                    imagecol_scaled = imagecol * expr@images[[i]]@scale.factors$lowres) %>%
      dplyr::inner_join(meta, by = "barcodeID")
    spatial_coord <- rbind(spatial_coord,tmp)
  }
  coords <- spatial_coord[,c("imagerow_scaled","imagecol_scaled")]
  colnames(coords) <- c('x','y')
  rownames(coords) <- spatial_coord$barcodeID
  coords <- coords[meta$barcodeID,]
  return(coords)
}

seurat.to.SpatialRNA <- function(Spatial) {
  nUMI = Spatial@meta.data$nCount_Spatial
  counts = Spatial[["Spatial"]]@counts
  names(nUMI) = colnames(counts)
  coords <- coords.file(Spatial)
  new("SpatialRNA", coords = coords, counts = counts, nUMI = nUMI)
}
seurat.to.scRNA <- function(sc.reference) {
  nUMI = sc.reference@meta.data$nCount_RNA
  counts = sc.reference[["RNA"]]@counts
  cell_types <-  Idents(sc.reference)
  names(nUMI) = colnames(counts)
  new("Reference", counts = counts,cell_types = cell_types,nUMI = nUMI)
}
######基因名称转换
pal2_ptr=read.delim('Pal_Ath_Ptr.txt',sep="\t",as.is=T,header=T,check.names=F)
ptr_genes=pal2_ptr[match(rownames(reference),pal2_ptr[,1]),]
rm=!duplicated(ptr_genes[,3])
reference_select=reference[rm,]
nUMI=reference_select@meta.data$nCount_RNA
counts = reference_select[["RNA"]]@counts
rownames(counts)=ptr_genes[rm,3]
cell_types <-  Idents(reference_select)
names(nUMI) = colnames(counts)
sc.reference=new("Reference", counts = counts,cell_types = cell_types,nUMI = nUMI)

# sc.reference <- seurat.to.scRNA(sc.reference)
SpatialRNA <-  seurat.to.SpatialRNA(spatial)
##NKT数量太多 用BHT201002_spacexr_fix.R替代
myRCTD <- create.RCTD(SpatialRNA, sc.reference, max_cores = 8,UMI_min = 1,CELL_MIN_INSTANCE = 25)
cell_type_names <- myRCTD@cell_type_info$info[[2]]
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
saveRDS(myRCTD,"myRCTD.RDS")
results <- myRCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
write.table(norm_weights,"cell_to_cluster_weights.xls",sep="\t")
write.table(results$results_df,"cell_to_cluster_type.xls",sep="\t")
Spatial <- readRDS('/data7/konghl/data/SpatialTranscriptome/BHT211041/analysis/combin.data.RDS')

metadata <- read.table("cell_to_cluster_weights.xls",row.names=1,header=T)
# colnames(metadata)=gsub('X','C',colnames(metadata))

# colnames(metadata)=
# metadata=metadata[,order(colnames(metadata))]  ### 细胞类型names排序
# colnames(metadata)=paste('C',colnames(metadata),sep="")
# colnames(metadata)=paste(colnames(metadata),sc_name,sep=":")

# Spatial=subset(Spatial,cells=rownames(metadata))
Spatial=spatial
for(i in colnames(metadata)){
  # Spatial[[i]]=0
  Spatial[[i]] = metadata[colnames(Spatial),i]
}

cell_types_all = colnames(metadata)
color <- c(getDefaultColors(n=length(cell_types_all)))
# samples = unique(Spatial[["orig.ident"]][,1])
samples=names(Spatial@images)
for (i in samples){
  print(i)
  p = scatterpie_plot1(Spatial,cell_types_all,slice =i,pie_scale=0.4)
  p1 <- p[[1]]
  bl <- p[[2]]
  # p2 = p1+scale_fill_manual(values=color)+guides(fill=FALSE)
  # ggsave(filename =  paste(i,"_heatmap_Nolegend.pdf",sep=""),p2, width = 10*bl, height = 10, dpi = 200)
  p3 = p1+scale_fill_manual(values=color,labels=)
  ggsave(filename =  paste(i,"_pieplot.pdf",sep=""),p3, width = 10*bl+1, height = 8.5, dpi = 200)
}
old=c('C1.0','C1.1','C1.2','C1.3','C1.4','C1.5','C11.0','C11.1','C14.0','C14.1')
new=c('C1.0 PIs','C1.1 PXy','C1.2 PPh','C1.3 PXy IN1','C1.4 MC',
      'C1.5 PPh IN1','C11.0 CZ','C11.1 DXy','C14.0 PLIs','C14.1 Ph')
