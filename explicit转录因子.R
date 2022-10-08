##########C12,C14,C11.0,C11.1,C8 get genecluster information
#######convert Populus geneID to arabidopsis TF gene
#######explicit analysis
### get genecluster path
library(dplyr)
library(tidyverse)
indir <- 'genecluster'
savePath <- getwd()
#### add gene annotation info
anno <- readRDS('anno_tf_go_kegg_ko.RDS')
#### read genecluster info
data_anno <- read.delim(file.path(indir,'branch_genecluster6_abcdef.txt'),sep="\t",as.is=T,header=T,check.names=F)
data_anno <- data_anno%>%tibble::column_to_rownames()%>%mutate(cluster=x,ID=anno[rownames(data_anno),'newGeneID'])
### data_anno=data.frame(gene=rownames(data_anno),cluster=data_anno[,1],ID=anno[rownames(data_anno),'newGeneID'])
# read arabidopsis TF gene
tf_info <- read.delim('Ath_TF_list.txt',sep="\t",as.is=T,header=T,check.names=F)
tf_id <- tf_info[,2]
########## filter Populus geneID
filter_gene_tf <- data_anno[data_anno[,'ID']%in%tf_id,]
filter_gene_tf <- cbind(filter_gene_tf,anno[filter_gene_tf[,'gene'],])
filter_gene_tf$gene <- as.character(filter_gene_tf$gene)
filter_gene_tf$ID <- as.character(filter_gene_tf$ID)

write.table(filter_gene_tf,file.path(savePath,'filter_tf.xlsx'),
            sep="\t",quote=F,row.names=F,col.names=T)
######### read seurat RDS 
expr <- readRDS('combin.subcluster.data.RDS')
Idents(expr) <- expr[['new_idents_labels']]

data_anno_tf_filter <- filter_gene_tf
rna <- expr[['Spatial']]@data
rna_filter <- rna[rownames(rna)%in%data_anno_tf_filter$gene,]
data_anno_tf_filter$new_id <- paste(data_anno_tf_filter[,'gene'],data_anno_tf_filter[,'cluster'],data_anno_tf_filter[,'ID'],sep="_")

rna1=rna_filter
rownames(data_anno_tf_filter)=data_anno_tf_filter[,'gene']
rownames(rna1)=data_anno_tf_filter[rownames(rna1),"new_id"]   

### target_data,tf_data,row is gene,column is samples
explicit<-function(target_data,tf_data,target_name,tf_name,path){
  library(MASS)
  filename=tf_name
  B=t(as.matrix(target_data)) # row samples col target_name
  # B=as.matrix(target_data)
  B=log2(B+0.00001)
  tf_data=t(as.matrix(tf_data))
  # tf_data=as.matrix(tf_data)
  tf_data_new=log2(tf_data+0.00001)
  A=cbind(1,tf_data_new)  #row samples col tf gene
  A=as.matrix(A)
  class(t(A))
  #[1] "matrix"
  x1=(t(A)%*%A)  #
  B=as.matrix(B)
  x2=t(A)%*%B
  # z=t(A)%*%B
  dim(x2)  ### row samples col targetgene
  #[1]  1679 29182
  beta=tryCatch({solve(x1)%*%x2},## If x1 is an invertible matrix，ginv(x1),MASS package
                error=function(x){library(MASS);ginv(x1)%*%x2})
  
  # library(MASS)
  # beta=ginv(x1)%*%x2
  
  Bp = A %*%beta  # row sample col target_gene
  c_sample=c()
  for(i in 1:nrow(Bp)){
    tmp=cor(Bp[i,],B[i,])
    if(i%%500==1){
      print(i)
    }
    c_sample=c(c_sample,tmp)
  }
  c=c_sample
  c_gene=c()
  for(i in 1:ncol(Bp)){
    tmp=cor(Bp[,i],B[,i])
    if(i%%500==1){
      print(i)
    }
    if(is.na(tmp)){
      print(i)
    }
    c_gene=c(c_gene,tmp)
  }
  cc=c_gene
  r=Bp-B
  rt=t(r)
  Bt=t(B)
  NRMSE=sqrt(sum(sum(rt^2))/sum(sum(Bt^2)))
  H=tryCatch({A%*%solve(x1)%*%t(A)},
             error=function(x){library(MASS);A%*%ginv(x1)%*%t(A)})
  
  ### unit matrix
  I = diag(nrow(A));
  J=matrix(rep(1,(nrow(Bp))^2),ncol=nrow(Bp))/nrow(Bp)
  SST=diag(t(B)%*%(I-J)%*%B)
  SSR=diag(t(B)%*%(H-J)%*%B)
  SSE=diag(t(B)%*%(I-H)%*%B)
  dof_sse = nrow(A) - ncol(A);
  dof_ssr = ncol(A) - 1;
  se_beta=tryCatch({sqrt(diag(solve(x1)) %*% t(SSE/dof_sse))},
                   error=function(x){sqrt(diag(ginv(x1)) %*% t(SSE/dof_sse))})
  tStat = beta/se_beta;
  beta_pvalue=pt(-abs(tStat),df=dof_sse)*2
  
  Fstat = (SSR / dof_ssr) / (SSE/dof_sse) 
  Fpvalue = pf(Fstat,dof_ssr,dof_sse);
  
  Fstat = Fstat
  qbeta = beta
  beta_pvalue = beta_pvalue
  SST = SST
  SSR = SSR
  SSE = SSE
  # idx = which(beta_pvalue <= 0.00001)
  idx=which(beta_pvalue<=1)
  e3 = beta[idx]
  e4 = beta_pvalue[idx];
  tf_name=c('intercept',tf_name)
  y=arrayInd(idx, c(nrow(beta),ncol(beta)))
  y1=y[,1]
  y2=y[,2]
  e1=tf_name[y[,1]]
  e2=target_name[y[,2]]
  SigEdges=data.frame(gene=e2,tf=e1,beta=e3,beta_pvalue=e4,stringsAsFactors = F)
  SigEdges=SigEdges[SigEdges[,2]!='intercept',]
  return(SigEdges)
}

out=c()
for(i in 1:nrow(rna1)){
  tmp_tf_name=rownames(rna1)[i]
  print(i)
  tmp_target_name=rownames(rna1)[-i]
  tmp_tf_data=t(data.frame(rna1[i,]))
  rownames(tmp_tf_data)=tmp_tf_name
  tmp_target_data=rna1[-i,]
  rownames(tmp_target_data)=tmp_target_name
  tmp_out=explicit(tmp_target_data,tmp_tf_data,tmp_target_name,tmp_tf_name,path)
  # out[[tmp_tf_name]]=tmp_out
  colnames(tmp_out)=c('target','tf','beta','beta_pvalue')
  out=rbind(out,tmp_out)
}
library(data.table)
out_copy=out
out[,c('target_gene','target_cluster','target_ID')]=do.call(cbind,tstrsplit(out[,'target'], "_", fixed=TRUE))
out[,c('tf_gene','tf_cluster','tf_ID')]=do.call(cbind,tstrsplit(out[,'tf'], "_", fixed=TRUE))

target_anno=anno[out$target_gene,]
colnames(target_anno)=paste('target',colnames(target_anno),sep="_")
tf_anno=anno[out$tf_gene,]
colnames(tf_anno)=paste('tf',colnames(tf_anno),sep="_")

results_pvalue_filter=results%>%filter(beta_pvalue<0.05)
results_pvalue_filter_beta_pvalue=results_pvalue_filter%>%filter(abs(beta)>0.3)
write.table(results_pvalue_filter,file.path(savePath,'results_pvalue_filter.xlsx'),sep="\t",quote=F,row.names=F,col.names=T)
write.table(results_pvalue_filter_beta_pvalue,file.path(savePath,'beta0.3_new.xlsx'),sep="\t",quote=F,row.names=F,col.names=T)
