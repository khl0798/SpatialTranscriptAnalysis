####
X2_mat=myRCTD@internal_vars_de$X2##X2是covariant值
#myRCTD@de_results$gene_fits$all_vals###是β*α的值
#X2*al_vals即为λ的值
#E[Yi| α] = λ¯i(α)eσ2ε,j /2
#而预测的期望值为λ*exp(σ^2)
#myRCTD@de_results$gene_fits$s_mat####standard errors for non-intercept terms
#myRCTD@de_results$gene_fits$con_mat#####check for covergence of each cell type
# means <- unlist(lapply(de_results_list, function(x) x$gene_fits$mean_val[gene,cell_type]))
### sds数据分布信息 gene_fits$I_mat 是方差？
#sds <- c(de_gene_results_08$gene_fits$I_mat[gene,ct_ind], de_gene_results_09$gene_fits$I_mat[gene,ct_ind], der_11$gene_fits$I_mat[gene,ct_ind])
#sds <- unlist(lapply(de_results_list, function(x) x$gene_fits$I_mat[gene,ct_ind]))
## myRCTD@de_results$gene_fits$mean_vals 是每个细胞类型对应的基因的logfc值


# Hello, to get cell type specific gene expression without differential expression, 
# you can run C-SIDE with only an intercept term. To get differential expression, 
# you should include at least one explanatory variable.
# After running DE, you can access the cell type specific gene expression values by looking at the 
# gene_fits$mean_val (only the DE parameter in case of single explanatory variable), 
# gene_fits$all_vals (all parameters) values, and the standard errors in gene_fits$s_mat.
#  To see the results for only the significant genes, look here: myRCTD@de_results$sig_gene_list. 
# More info here https://github.com/dmcable/spacexr/tree/master/documentation.


cell_type <- 'CAF'
log_fc <- myRCTD@de_results$gene_fits$mean_val[gene_list_type,cell_type]  ### mean_val对应的基因的logfc值
cell_type_ind <- which(myRCTD@internal_vars_de$cell_types == cell_type)*2  ###
z_score <- log_fc / myRCTD@de_results$gene_fits$I_mat[gene_list_type, cell_type_ind] ####I_mat对应的方差
### z_score的值为log_fc/I_mat
p_vals <- 2*(1-pnorm(abs(z_score))) ##计算p值
p_vals[p_vals == 0] <- 1e-16

plot_df <- data.frame(log_fc, -log(p_vals,10))
colnames(plot_df) <- c('log_fc','p_val')
MAX_P_VAL <- -log(max(res_genes$p_val),10)
plot_df$gene <- rownames(plot_df)
#plot_df$label <- plot_df$gene
#plot_df$label[!(plot_df$label %in% rownames(res_genes))] <- NA
plot_df$sig <- (rownames(plot_df) %in% rownames(res_genes))
plot_df$log_fc <- plot_df$log_fc*log(exp(1),2)    ###转成log2fc的值
plot_df$group <- 'reg'

spacexr:::predict_CSIDE<-function (cell_type_ind, gene_fits, gene, X2_mat){
    sigma <- as.numeric(gene_fits$sigma_g[gene])/100
    predictions <- exp(X2_mat %*% gene_fits$all_vals[gene, , cell_type_ind])
    return(predictions * exp(sigma^2/2))
}
plot_prediction_genes <- function(cell_type, barcodes, my_beta, puck, X2, sig_genes, doublet_mode, datadir, gene_fits) {
  if(!dir.exists(file.path(datadir,'prediction_plots')))
    dir.create(file.path(datadir,'prediction_plots'))
  MULT <- 500
  cell_type_ind <- which(cell_type == colnames(my_beta))
  gene_list_sig <- rownames(sig_genes)
  if(!doublet_mode)
    sing_thresh <- 0.8
  else
    sing_thresh <- 0.999
  barcodes_sing <- names(which(my_beta[barcodes,cell_type] > sing_thresh))
  plots <- list()
  if(length(gene_list_sig) > 0 & length(barcodes_sing) > 0) {
    for(i in 1:length(gene_list_sig)) {#length(gene_list_sig)
      gene = gene_list_sig[i]
      predictions <- predict_GLAMDE(cell_type_ind, gene_fits, gene, X2[barcodes_sing,])[,1]
      plots[[i]] <- plot_puck_continuous(puck, barcodes_sing, predictions*MULT,
                                         ylimit = c(0,quantile(predictions*MULT, 0.95)), title = gene)
    }
    pdf(file.path(datadir,paste0("prediction_plots/de_genes_",stringr::str_replace(cell_type,'/','_'),".pdf")))
    invisible(lapply(plots, print))
    dev.off()
  }
}
###画图基因的表达期望值
#' Makes a spatial plot of GLAMDE fitted gene expression
#'
#' Units counts per 500
#'
#' @param myRCTD \code{\linkS4class{RCTD}} object containing \code{de_results}, after running GLAMDE
#' @param gene gene to be plotted
#' @param cell_type cell_type to be plotted (only single cell type pixels)
#' @return plot of fitted gene expression
#' @export
plot_prediction_gene <- function(myRCTD, cell_type, gene) {
  barcodes <- myRCTD@internal_vars_de$barcodes ##挑选的barcodes
  my_beta <- myRCTD@internal_vars_de$my_beta##即weights值
  puck <- myRCTD@spatialRNA 
  doublet_mode <- myRCTD@internal_vars_de$doublet_mode
  cell_type_ind <- which(colnames(my_beta) == cell_type)
  MULT <- 500
  if(!doublet_mode)
    sing_thresh <- 0.8
  else
    sing_thresh <- 0.999
  barcodes_sing <- names(which(my_beta[barcodes,cell_type] > sing_thresh))
  plots <- list()
  if(length(barcodes_sing) >= 10) {
    predictions <- predict_GLAMDE(cell_type_ind, myRCTD@de_results$gene_fits, gene,
                                 myRCTD@internal_vars_de$X2[barcodes_sing,])[,1]
    p <- plot_puck_continuous(puck, barcodes_sing, predictions*MULT,
                                       ylimit = c(0,quantile(predictions*MULT, 0.95)), title = gene)
    return(p)
  } else {
    print('plot_prediction_gene: not plotting because fewer than 10 cell type singlet pixels found')
  }
}
#####故意设置差异表达基因
#change_gene <- 'Lsamp'
change_genes=rownames(myRCTD@spatialRNA@counts)[1:50]
high_barc <- names(explanatory.variable[explanatory.variable > 0.5])
low_barc <- names(explanatory.variable[explanatory.variable < 0.5])
for(change_gene in change_genes){
myRCTD@originalSpatialRNA@counts[change_gene, high_barc] <- myRCTD@spatialRNA@counts[change_gene, high_barc] * 3 
}
myRCTD@config$max_cores <- 20
myRCTD <- run.CSIDE.single(myRCTD, explanatory.variable, gene_threshold = .00001, 
                        cell_type_threshold = 125, fdr =0.1,doublet_mode = F) 

###########重新根据
predict_GLAMDE<-function (cell_type_ind, gene_fits, gene, X2_mat){
  sigma <- as.numeric(gene_fits$sigma_g[gene])/100
  predictions <- exp(X2_mat %*% gene_fits$all_vals[gene, , cell_type_ind])
  return(predictions * exp(sigma^2/2))
}
cur_cell_types <- myRCTD@internal_vars_de$cell_types
x=get_quant_df(myRCTD, myRCTD@de_results$gene_fits, cell_types, cur_cell_types, gene="IGKC", 
               multi_region = F, prop_thresh = 0.999, param_position = 2)

get_quant_df<-function (RCTDde, gene_fits, cell_types, cur_cell_types, gene, 
                        multi_region = F, prop_thresh = 0.999, param_position = 2) 
{
  my_beta <- RCTDde@internal_vars_de$my_beta
  pred_tot <- numeric(dim(my_beta)[1])
  X2 <- RCTDde@internal_vars_de$X2[rownames(my_beta), param_position]
  if (multi_region) 
    X2 <- apply(RCTDde@internal_vars_de$X2, 1, function(x) which(as.logical(x)))
  nUMI <- RCTDde@spatialRNA@nUMI[rownames(my_beta)]
  Y <- RCTDde@originalSpatialRNA@counts[gene, rownames(my_beta)]/nUMI  ###线性归一化的结果
  for (cell_type_ind in 1:length(cell_types)) {
    sigma <- as.numeric(gene_fits$sigma_g[gene])/100
    pred_ct <- predict_GLAMDE(cell_type_ind, gene_fits, gene, 
                              RCTDde@internal_vars_de$X2[rownames(my_beta), ])
    pred_tot <- pred_tot + pred_ct * my_beta[, cell_type_ind]
  }
  a_pred <- pred_tot/exp(sigma^2/2)
  var_pred <- a_pred^2 * exp(sigma^2/2) * (exp(sigma^2/2) - 
                                             1)
  var_tot <- var_pred + pred_tot/nUMI
  if (length(cur_cell_types) > 1) 
    my_barc <- which(rowSums(my_beta[, cur_cell_types]) > 
                       prop_thresh)
  else my_barc <- which(my_beta[, cur_cell_types] > prop_thresh)
  if (!multi_region) {
    plot_df <- cbind(pred_tot[my_barc], var_tot[my_barc], 
                     my_beta[my_barc, min(which(cell_types %in% cur_cell_types))], 
                     X2[my_barc], Y[my_barc])
    colnames(plot_df) <- c("pred", "var", "berg_weight", 
                           "region", "Y")
  }else {
    plot_df <- cbind(pred_tot[my_barc], var_tot[my_barc], 
                     my_beta[my_barc, cur_cell_types], X2[my_barc], Y[my_barc])
    colnames(plot_df) <- c("pred", "var", cur_cell_types, 
                           "region", "Y")
  }
  plot_df <- data.frame(plot_df)
}

get_de_pop <- function(cell_type, de_results_list) {
  ct_ind <- which(colnames(de_results_list[[1]]$gene_fits$mean_val) == cell_type)*2
  #gene_list <- intersect(intersect(rownames(de_gene_results_08$gene_fits$mean_val), rownames(de_gene_results_09$gene_fits$mean_val)),rownames(der_11$gene_fits$mean_val))
  gene_list <- Reduce(intersect, lapply(de_results_list, function(x) names(which(x$gene_fits$con_mat[,cell_type]))))
  de_pop <- matrix(0, nrow = length(gene_list), ncol = 4)
  colnames(de_pop) <- c('sig_p', 'mean_est', 'sd_est', 'Z_est')
  rownames(de_pop) <- gene_list
  for(gene in gene_list) {
    #means <- c(de_gene_results_08$gene_fits$mean_val[gene,cell_type], de_gene_results_09$gene_fits$mean_val[gene,cell_type], der_11$gene_fits$mean_val[gene,cell_type])
    means <- unlist(lapply(de_results_list, function(x) x$gene_fits$mean_val[gene,cell_type]))
    #sds <- c(de_gene_results_08$gene_fits$I_mat[gene,ct_ind], de_gene_results_09$gene_fits$I_mat[gene,ct_ind], der_11$gene_fits$I_mat[gene,ct_ind])
    sds <- unlist(lapply(de_results_list, function(x) x$gene_fits$I_mat[gene,ct_ind]))
    sds[is.na(sds)] <- 1000
    sig_p <- sqrt(estimate_var(means, sds))
    var_t <- sds^2 + sig_p^2
    var_est <- 1/sum(1 / var_t)
    mean_est <- sum(means / var_t)*var_est
    
    sd_est <- sqrt(var_est)
    Z_est <- mean_est / sd_est
    de_pop[gene, ] <- c(sig_p, mean_est, sd_est, Z_est)
  }
  de_pop <- as.data.frame(de_pop)
  return(de_pop)
}

# one_ct_genes <- function(cell_type, myRCTD_list, de_results_list, resultsdir, cell_types_present, q_thresh = .01, p_thresh = 1, filter = T, order_gene = F, plot_results = T) {
#   de_pop <- get_de_pop(cell_type, de_results_list)
#   myRCTD <- myRCTD_list[[1]]
#   gene_big <- Reduce(intersect, lapply(myRCTD_list, 
#                               function(myRCTD) get_gene_list_type_wrapper(myRCTD, 
#                               cell_type, cell_types_present)))
#   cell_type_means <- myRCTD@cell_type_info$info[[1]][gene_big,cell_types_present]
#   cell_prop <- sweep(cell_type_means,1,apply(cell_type_means,1,max),'/')
#   p_vals <- 2*(1-pnorm(abs(de_pop[gene_big,'Z_est'])))
#   names(p_vals) <- gene_big
#   q_vals<- p.adjust(p_vals,'BH')
#   if(filter)
#     gene_final <- intersect(gene_big[which(q_vals < q_thresh & p_vals < p_thresh)],
#                        gene_big[which(abs(de_pop[gene_big,'mean_est']) > 0.4)])
#   else
#     gene_final <- gene_big
#   final_df <- cbind(de_pop[gene_final,],cell_prop[gene_final,c(cell_type)], 
#                     cell_type_means[gene_final,cell_type], q_vals[gene_final])
#   colnames(final_df) <- c( 'sig_p',    'mean_est'  ,   'sd_est'    ,  'Z_est'    , 'ct_prop' ,'expr' ,'q_val')
#   final_df$p <- 2*(1 - pnorm(abs(final_df$Z_est)))
#   L <- length(myRCTD_list)
#   mean_sd_df <- matrix(0, nrow = length(gene_final), ncol = L*2)
#   rownames(mean_sd_df) <- gene_final
#   colnames(mean_sd_df) <- c(unlist(lapply(1:L, function(x) paste('mean', x))), unlist(lapply(1:L, function(x) paste('sd', x))))
#   for(gene in gene_final) {
#     m_sd <- get_means_sds(cell_type, gene, de_results_list)
#     mean_sd_df[gene,] <- c(m_sd$means, m_sd$sds)
#   }
#   final_df <- cbind(final_df, mean_sd_df)
#   if(order_gene)
#     final_df <- final_df[order(gene_final), ]
#   else
#     final_df <- final_df[order(-abs(final_df$mean_est)),]
#   #plot(log(final_df$expr,10), log(final_df$p,10))
#   if(plot_results) {
#     print('writing')
#     write.csv(final_df,file.path(resultsdir,paste0(cell_type,'_cell_type_genes.csv')))
#   }
#   print('done')
#   return(list(de_pop = de_pop, gene_final = gene_final))
# }
spacexr:::one_ct_genes<-function (cell_type, myRCTD_list, de_results_list, resultsdir, 
    cell_types_present, q_thresh = 0.01, p_thresh = 1, filter = T, 
    order_gene = F, plot_results = T, use.groups = F, group_ids = NULL, 
    MIN.CONV.REPLICATES = 2, MIN.CONV.GROUPS = 2, CT.PROP = 0.5, 
    log_fc_thresh = 0.4, normalize_expr = F) 
{
    print(paste0("one_ct_genes: population inference on cell type, ", 
        cell_type))
    myRCTD <- myRCTD_list[[1]]
    cell_type_means <- myRCTD@cell_type_info$info[[1]][, cell_types_present]
    cell_prop <- sweep(cell_type_means, 1, apply(cell_type_means, 
        1, max), "/")
    de_pop <- get_de_pop(cell_type, de_results_list, cell_prop, 
        use.groups = use.groups, group_ids = group_ids, MIN.CONV.REPLICATES = MIN.CONV.REPLICATES, 
        MIN.CONV.GROUPS = MIN.CONV.GROUPS, CT.PROP = CT.PROP)
    gene_big <- rownames(de_pop)[which(de_pop$tau >= 0)]
    p_vals <- 2 * (pnorm(-abs(de_pop[gene_big, "Z_est"])))
    names(p_vals) <- gene_big
    q_vals <- p.adjust(p_vals, "BH")
    if (filter) 
        gene_final <- intersect(gene_big[which(q_vals < q_thresh & 
            p_vals < p_thresh)], gene_big[which(abs(de_pop[gene_big, 
            "log_fc_est"]) > log_fc_thresh)])
    else gene_final <- gene_big
    gene_df <- cbind(de_pop[gene_big, ], cell_prop[gene_big, 
        c(cell_type)], cell_type_means[gene_big, cell_type], 
        q_vals[gene_big])
    colnames(gene_df) <- c(colnames(de_pop), "ct_prop", "expr", 
        "q_val")
    gene_df$p <- 2 * (pnorm(-abs(gene_df$Z_est)))
    final_df <- gene_df[gene_final, ]
    L <- length(myRCTD_list)
    mean_sd_df <- matrix(0, nrow = length(gene_final), ncol = L * 
        2)
    rownames(mean_sd_df) <- gene_final
    colnames(mean_sd_df) <- c(unlist(lapply(1:L, function(x) paste("mean", 
        x))), unlist(lapply(1:L, function(x) paste("sd", x))))
    for (gene in gene_final) {
        m_sd <- get_means_sds(cell_type, gene, de_results_list)
        mean_sd_df[gene, ] <- c(m_sd$means, m_sd$sds)
    }
    final_df <- cbind(final_df, mean_sd_df)
    if (length(gene_final) > 1) 
        if (order_gene) 
            final_df <- final_df[order(gene_final), ]
        else final_df <- final_df[order(-abs(final_df$log_fc_est)), 
            ]
    if (plot_results) {
        print("writing")
        write.csv(final_df, file.path(resultsdir, paste0(cell_type, 
            "_cell_type_genes.csv")))
    }
    print("done")
    return(list(de_pop = gene_df, gene_final = gene_final, final_df = final_df))
}
spacexr:::get_de_pop<-function (cell_type, de_results_list, cell_prop, use.groups = F, 
    group_ids = NULL, MIN.CONV.REPLICATES = 2, MIN.CONV.GROUPS = 2, 
    CT.PROP = 0.5, S.MAX = 4) 
{
    if (!use.groups) 
        group_ids <- NULL
    ct_ind <- which(colnames(de_results_list[[1]]$gene_fits$mean_val) == 
        cell_type) * 2
    gene_list <- Reduce(union, lapply(de_results_list, function(x) names(which(x$gene_fits$con_mat[, 
        cell_type]))))
    gene_list <- intersect(gene_list, rownames(cell_prop)[(which(cell_prop[, 
        cell_type] >= CT.PROP))])
    if (!use.groups) {
        de_pop <- matrix(0, nrow = length(gene_list), ncol = 5)
        colnames(de_pop) <- c("tau", "log_fc_est", "sd_est", 
            "Z_est", "p_cross")
    }
    else {
        group_names <- unique(group_ids)
        n_groups <- length(group_names)
        de_pop <- matrix(0, nrow = length(gene_list), ncol = 6 + 
            2 * n_groups)
        colnames(de_pop) <- c("tau", "log_fc_est", "sd_est", 
            "Z_est", "p_cross", "delta", unlist(lapply(group_names, 
                function(x) paste0(x, "_group_mean"))), unlist(lapply(group_names, 
                function(x) paste0(x, "_group_sd"))))
    }
    rownames(de_pop) <- gene_list
    ii <- 1
    for (gene in gene_list) {
        ii <- ii + 1
        if (ii%%1000 == 0) 
            message(paste("get_de_pop: testing gene,", gene, 
                ", of index:", ii))
        check_con <- function(x) {
            ifelse(gene %in% rownames(x$gene_fits$con_mat), x$gene_fits$con_mat[gene, 
                cell_type] && !is.na(x$gene_fits$s_mat[gene, 
                ct_ind]) && (x$gene_fits$s_mat[gene, ct_ind] < 
                S.MAX), FALSE)
        }
        con <- unlist(lapply(de_results_list, check_con))
        if (use.groups) 
            con <- unname(con & table(group_ids[con])[as.character(group_ids)] >= 
                2)
        used_groups <- names(table(group_ids[con]))
        if (sum(con) < MIN.CONV.REPLICATES || (use.groups && 
            length(used_groups) < MIN.CONV.GROUPS)) {
            if (use.groups) 
                de_pop[gene, ] <- c(-1, 0, 0, 0, 0, 0, rep(0, 
                  n_groups), rep(-1, n_groups))
            else de_pop[gene, ] <- c(-1, 0, 0, 0, 0)
        }
        else {
            means <- unlist(lapply(de_results_list[con], function(x) x$gene_fits$mean_val_cor[[cell_type]][gene]))
            sds <- unlist(lapply(de_results_list[con], function(x) x$gene_fits$s_mat[gene, 
                ct_ind]))
            sds[is.na(sds)] <- 1000
            if (is.null(group_ids)) 
                gid <- NULL
            else gid <- group_ids[con]
            sig_p <- estimate_tau_group(means, sds, group_ids = gid)
            var_t <- sds^2 + sig_p^2
            if (!use.groups) {
                var_est <- 1/sum(1/var_t)
                mean_est <- sum(means/var_t) * var_est
                p_cross <- get_p_qf(means, sds)
            }
            else {
                S2 <- 1/(aggregate(1/var_t, list(group_ids[con]), 
                  sum)$x)
                E <- (aggregate(means/var_t, list(group_ids[con]), 
                  sum)$x) * S2
                Delta <- estimate_tau_group(E, sqrt(S2))
                var_T <- (Delta^2 + S2)
                var_est <- 1/sum(1/var_T)
                mean_est <- sum(E/var_T) * var_est
                p_cross <- get_p_qf(E, sqrt(S2))
                E_all <- rep(0, n_groups)
                s_all <- rep(-1, n_groups)
                names(E_all) <- group_names
                names(s_all) <- group_names
                E_all[used_groups] <- E
                s_all[used_groups] <- sqrt(S2)
            }
            sd_est <- sqrt(var_est)
            Z_est <- mean_est/sd_est
            if (use.groups) 
                de_pop[gene, ] <- c(sig_p, mean_est, sd_est, 
                  Z_est, p_cross, Delta, E_all, s_all)
            else de_pop[gene, ] <- c(sig_p, mean_est, sd_est, 
                Z_est, p_cross)
        }
    }
    de_pop <- as.data.frame(de_pop)
    return(de_pop)
}
estimate_var <- function(x, s) {
  return(max(var(x) - mean(s^2),0))
}

########解析表达值
barcodes <- intersect(names(myRCTD@spatialRNA@nUMI), colnames(myRCTD@spatialRNA@counts))
X <- as.matrix(rep(1, length(barcodes)))
rownames(X) <- barcodes
myRCTD_expr <- run.CSIDE(myRCTD, X, barcodes,params_to_test=1)

#############
spacexr:::predict_CSIDE_all<-function (RCTDde, gene){
    cell_types <- RCTDde@internal_vars_de$cell_types
    gene_fits <- RCTDde@de_results$gene_fits
    my_beta <- RCTDde@internal_vars_de$my_beta
    pred_tot <- numeric(dim(my_beta)[1])
    sigma <- as.numeric(gene_fits$sigma_g[gene])/100
    for (cell_type_ind in 1:length(cell_types)) {
        pred_ct <- predict_CSIDE(cell_type_ind, gene_fits, gene, 
            RCTDde@internal_vars_de$X2[rownames(my_beta), ])
        pred_tot <- pred_tot + pred_ct * my_beta[, cell_type_ind]
    }
    return(pred_tot)
}



