# SNF patient similarity



# Settings
suppressMessages(library(dplyr))
suppressMessages(library(SNFtool))
# suppressMessages(library(CancerSubtypes))

suppressMessages(library(Mfuzz))
suppressMessages(library(e1071))




#' Custimize SNF analysis pipeline
#'
#' @param dat list of all kinds of data matrix
#' @param K 
#' @param alpha 
#' @param dist_metric 
#' @param min_cluster_num 
#' @param max_cluster_num 
#' @param cluster_method 
#' @param display_cluster 
#' @param workdir 
#' @param prefix 
#'
#' @return
#' @export
#'
#' @examples
do_snf_analysis = function(dat=list(),
                           K=20,
                           alpha=0.5,
                           T=10,
                           dist_metric="euclidean",
                           min_cluster_num=2,
                           max_cluster_num=10,
                           cluster_method="spectral",
                           display_cluster=FALSE,
                           workdir="./",
                           prefix="snf_graph"){
  
  # K : number of neighbors, usually (10~30)
  # alpha : hyperparameter, usually (0.3~0.8)
  # T : Number of Iterations, usually (10~20)
  # C : number of clusters
  ## Input data is of size n x d_1, where n is the number of patients, d_1 is the number of genes.
  
  print("Perform SNF clustering analysis")
  ## If the data are all continuous values, we recommend the users to perform standard normalization before using SNF, 
  ## though it is optional depending on the data the users want to use.
  message("*** do standard normalization")
  # use SNF with multiple views
  dat_norm = lapply(dat, function(x) SNFtool::standardNormalization(x)) 
  
  ## Calculate distance matrices (here we calculate Euclidean Distance, you can use other distance, e.g,correlation)
  ## Calculate the pair-wise distance; If the data is continuous, we recommend to use the function "dist2" as follows; 
  # if the data is discrete, we recommend the users to use ""
  message("*** calculate distance matrices")
  if(dist_metric=="euclidean"){
    dat_dist = lapply(dat_norm, function(x) SNFtool::dist2(as.matrix(x), as.matrix(x))^(1/2))
  }else{
    dat_dist = lapply(dat_norm, function(x) cor(as.matrix(x), as.matrix(x), method = "pearson")) # correlation
  }
  
  ## next, construct similarity graphs 
  message("*** construct similarity graphs")
  dat_graph = lapply(dat_dist, function(x) affinityMatrix(x, as.numeric(K), as.numeric(alpha)))
  
  ## next, fuse all the graphs
  ## then the overall matrix can be computed by similarity network fusion(SNF):
  message("*** fuse all the graphs")
  fusion_graph = SNF(dat_graph, K = K, t = T) # t: Number of iterations for the diffusion process.
  write.table(fusion_graph, file=file.path(workdir,paste0(prefix,"_fusion_graph.tsv")), row.names = TRUE, quote=F, sep = "\t")
  
  dat_graph[["fusion_graph"]] = fusion_graph
  
  ## With this unified graph W (fusion_graph) of size n x n, you can do either spectral clustering or Kernel NMF. 
  ## working on Spectral Clustering
  message("*** work on Spectral Clustering")
  smp_lst_fused_lst = list()
  for(C in seq(min_cluster_num,max_cluster_num)){
    message(paste0("      clustering samples into ",C," groups"))
    if(cluster_method=="spectral"){
      group = spectralClustering(fusion_graph, C) 	# the final subtypes information
    }else if(cluster_method=="kmeans"){
      # m1 = Mfuzz::mestimate(as.data.frame(fusion_graph))
      # group = Mfuzz::mfuzz(fusion_graph, c = C, m=1.25)
      
    }else if(tolower(cluster_method)=="nmf"){
      kernel_result = do_kernel_nmf_clustering(matrix_x = fusion_graph, 
                                               matrix_y = NULL, 
                                               sigma = 0.05, 
                                               C = C, 
                                               plot = FALSE, 
                                               workdir = "./", 
                                               prefix = "knmf") 
    }
    
    cluster_result = data.frame(cluster=group, 
                                Case_ID=colnames(fusion_graph))
    write.table(cluster_result, file.path(workdir,paste0(prefix,"_",C,"cluster_result.tsv")),row.names = F,sep="\t",quote=F)
    
    if(display_cluster){
      png(file=file.path(workdir,paste0(prefix,"_",C,"cluster.png")), height = 8, width = 8, units = "in", res = 300)
      displayClustersWithHeatmap(fusion_graph, group)
      dev.off()
    }

    smp_lst_fused_lst[[paste0(C,"_clusters")]] = group
    
  }
  
  return(list(dat_norm=dat_norm,
              dat_dist=dat_dist,
              dat_graph=dat_graph,
              group=smp_lst_fused_lst))
}



#' Evaluate SNF clustering results
#'
#' @param fusion_graph 
#' @param group_list 
#' @param dat_dist 
#' @param plot.sil 
#' @param plot.pca 
#' @param plot.heatmap 
#' @param meta_df 
#' @param workdir 
#' @param prefix 
#'
#' @return
#' @export
#'
#' @examples
evaluate_snf_cluster = function(fusion_graph,
                                group_list=NULL,
                                dat_list=NULL,
                                plot.pca=FALSE,
                                plot.umap=FALSE,
                                plot.sil=FALSE,
                                plot.heatmap=FALSE,
                                binary.ht=FALSE,
                                meta_df=NULL,
                                arm_annot=NULL,
                                mut_col=NULL,
                                custom_color=FALSE,
                                myColor=c("red","blue"),
                                workdir="./",
                                prefix="snf_cluster"){
  
  print("Evaluate SNF clustering results")
  ## With this unified graph W of size n x n, you can do either spectral clustering or Kernel NMF. 
  ## for example, spectral clustering
  source("../cluster/AML_heatmap_column_annotation.R")
  source("../cluster/silhouette_analysis.R")

  # color_list = c("tomato", "forest green", "dark blue", "purple2", "goldenrod4", "gray20","darkcyan","darkred","#FF9083","#F90783")
  
  for(C in gsub("_clusters","",names(group_list))){
    C = as.numeric(C)
    group = group_list[[paste0(C,"_clusters")]]
    
    # get sample names in clustering result, implement from displayClustersWithHeatmap
    W = fusion_graph
    normalize <- function(X) X/rowSums(X)
    ind <- sort(as.vector(group), index.return = TRUE)
    sorted_group = ind$x
    ind <- ind$ix

    diag(W) <- median(as.vector(W))
    
    W <- normalize(W)
    W <- W + t(W)
    W_ht = as.data.frame(W)[ind, ind]
    # colnames(W_ht) = sapply(colnames(W_ht), function(x) gsub("[.]","-",x))
    case_ids = colnames(W_ht)
    
    
    group_list[[paste0(C,"_clusters")]] = data.frame(Cluster=sorted_group, #group_list[[paste0(C,"_clusters")]],
                                                     Case_ID=case_ids)
                                                     
    
    if(plot.pca){
      print(paste0("*** plot PCA for ", C))
      pca_plot = plot_pca(W = fusion_graph, group = group, C = C)
      ggsave(filename=paste0(file.path(workdir,paste0("pca_",prefix,"_",C,"cluster.png"))),
             plot=pca_plot,
             width = 8,
             height = 8,
             dpi = 300,
             units = "in")
    }
    
    if(plot.umap){
      print(paste0("*** plot UMAP for ", C))
      umap_plot = plot_umap(W = fusion_graph, 
                            C = C, 
                            group = group, 
                            add_ellipse = FALSE, 
                            n_components = 2, 
                            n_neighbors = 20,
                            add_segments = TRUE,
                            workdir = workdir, 
                            prefix = C)
      
      umap_plot
      ggsave(filename=paste0(file.path(workdir,paste0("umap_",prefix,"_",C,"cluster.png"))),
             # plot=umap_plot,
             width = 8,
             height = 8,
             dpi = 300,
             units = "in")

    }
    
    if(plot.sil && plot.umap){
      print(paste0("*** plot silhouette score for ", C))
      umap_result = calculate_umap_dist(filename = file.path(workdir, paste0(C,"_umap.tsv")))
      umap_dist = umap_result[["dist"]]
      umap_group = umap_result[["group"]]
      sil_result = plot_silhouette(silhouette(umap_group, umap_dist))
      sil_plot = sil_result[["plot"]]
      sil_df = sil_result[["sil_df"]]
      write.table(sil_df, file = file.path(workdir,paste0("sil_snf_",C,"clusters.tsv")),sep = "t",quote = F,row.names = F)
      
      # sil_plot
      ggsave(sil_plot,
             filename=file.path(workdir,paste0("sil_snf_",C,"clusters.png")),
             width = 8,
             height = 8,
             dpi = 300,
             units = "in")
    }
    
    if(plot.heatmap){
      print(paste0("*** plot Heatmap for ", C))
      # prepare annotation matrix for heatmap tracks
      # column_ha = aml_ht_column_annotation(meta_df, 
      #                                      case_ids,
      #                                      case_col = 'Case_ID')
      # print(head(case_ids))
      column_ha = aml_ht_column_annotation.2alec(meta_df,
                                                 arm=arm_annot,
                                                 case_col="patientID",
                                                 samples=case_ids,
                                                 # cluster_idx=ind$x,
                                                 mut_col=mut_col,
                                                 na_col="lightgrey")

      if(binary.ht){
        mat = matrix(0, nrow = nrow(W_ht), ncol = ncol(W_ht))
        rownames(mat) = rownames(W_ht)
        colnames(mat) = colnames(W_ht)
        for(i in sort(unique(sorted_group), decreasing = F)){
          idx = which(sorted_group==i)
          idx_start = idx[1]
          idx_end = idx[length(idx)]
          if(custom_color){
            mat[idx_start:idx_end,idx_start:idx_end] = i
          }else{
            mat[idx_start:idx_end,idx_start:idx_end] = 1
          }
        }
        if(custom_color){
          paletteLength = C+1
          myColor = viridis::viridis(paletteLength)
          myColor = c("grey",myColor)
        }
        
        ht1 = suppressMessages(Heatmap(mat, 
                                       name = paste0(C,"_clusters"),
                                       # show_column_dend = TRUE,
                                       # show_row_dend = T,
                                       #row_labels = c(sel_prot_mito$ID, rep("", 3000)),
                                       cluster_rows = FALSE,
                                       cluster_columns = FALSE,
                                       # cluster_columns = cluster_within_group(mat, sorted_group),
                                       # row_split = sorted_group,
                                       column_split = sorted_group,
                                       top_annotation = column_ha, 
                                       show_row_names = FALSE,
                                       col = myColor,
                                       # row_split = 3,
                                       # heatmap_legend_param = list(color_bar="continuous", #at=c(-2,-1,0,1,2),
                                       #                             legend_direction="vertical", 
                                       #                             legend_height=unit(4,"cm"),
                                       #                             title_position="topcenter", 
                                       #                             title_gp=gpar(fontsize=10, fontface="bold")),
                                       column_names_gp = grid::gpar(fontsize = 4),
                                       row_names_gp = grid::gpar(fontsize = 4)))

      }else{
        mat = as.matrix(W_ht)
        
        ht1 = suppressMessages(Heatmap(mat, 
                                       name = paste0(C,"_clusters"),
                                       # show_column_dend = TRUE,
                                       # show_row_dend = T,
                                       #row_labels = c(sel_prot_mito$ID, rep("", 3000)),
                                       cluster_rows = FALSE,
                                       cluster_columns = FALSE,
                                       # cluster_columns = cluster_within_group(mat, sorted_group),
                                       # row_split = sorted_group,
                                       column_split = sorted_group,
                                       top_annotation = column_ha, 
                                       show_row_names = FALSE,
                                       # row_split = 3,
                                       # heatmap_legend_param = list(color_bar="continuous", #at=c(-2,-1,0,1,2),
                                       #                             legend_direction="vertical", 
                                       #                             legend_height=unit(4,"cm"),
                                       #                             title_position="topcenter", 
                                       #                             title_gp=gpar(fontsize=10, fontface="bold")),
                                       column_names_gp = grid::gpar(fontsize = 4),
                                       row_names_gp = grid::gpar(fontsize = 4)))
      }

      png(file.path(workdir,paste0("heatmap_",prefix,"_",C,"cluster.png")), width=15, height=15, units="in", res=300) 
      draw(ht1, merge_legends = TRUE, 
           heatmap_legend_side="left", 
           annotation_legend_side="left",
           legend_grouping = "original")
      
      dev.off()
    }
  }

  return(group_list)
}



#' Get graph concordance by using NMI to measure the goodness of the obtained labels. 
#'
#' @param fusion_graph 
#' @param dat_graph 
#' @param sample_annot 
#' @param cluster_num 
#' @param workdir 
#' @param prefix 
#' @param plot 
#'
#' @return
#' @export
#'
#' @examples
evaluate_graph_concordance = function(fusion_graph,
                                      dat_graph,
                                      sample_annot=NULL,
                                      min_cluster_num=2,
                                      max_cluster_num=10,
                                      workdir="./",
                                      prefix="snf_cluster",
                                      plot=FALSE){
  # ## you can evaluate the goodness of the obtained clustering results by calculate Normalized mutual information (NMI): 
  # ## if NMI is close to 1, it indicates that the obtained clustering is very close to the "true" cluster information; 
  # ## if NMI is close to 0, it indicates the obtained clustering is not similar to the "true" cluster information.
  if(!is.null(sample_annot)){
    png(file=file.path(workdir,"snf_cluster_goodness.png"), width = 8, height = 8, res=300, units = "in")
    displayClusters(fusion_graph, group)
    dev.off()
    
    SNFNMI = calNMI(group, sample_annot)
    write.table(SNFNMI, file=file.path(workdir,"snf_NMI.tsv"), row.names = F, quote=F, sep="\t")
  }
  
  ## you can also find the concordance between each individual network and the fused network
  # ConcordanceMatrix = concordanceNetworkNMI(c(list(Fused_all=fusion_graph), dat_graph), C=cluster_num)
  # Note: Concordance is evaluated based on normalized mutual information (NMI)
  # The similarity between two community structures is quantified by the Normalized Mutual Information (NMI) (Alexander-Bloch et al. 2012; Kuncheva and Hadjitodorov 2004; Strehl and Ghosh 2003), which is one of the most popular, similarity measure based on Taya et al. Applied Network Science (2016) 1:8 Page 7 of 2
  # information theory.
  message(paste0("*** find the concordance between each individual network and the fused network, concordance or similarity is quantified by the Normalized Mutual Information (NMI)"))
  
  nmi_list = list()
  for(cluster_num in seq(min_cluster_num, max_cluster_num)){
    ConcordanceMatrix = suppressWarnings(concordanceNetworkNMI(c(list(Fusion=fusion_graph), dat_graph), C=cluster_num))
    rownames(ConcordanceMatrix) = names(c(list(Fusion=fusion_graph), dat_graph))
    colnames(ConcordanceMatrix) = names(c(list(Fusion=fusion_graph), dat_graph))
    message(paste0("       working on ",cluster_num," clusters setting"))
    
    # print(ConcordanceMatrix)
    if(plot){
      source("../shared/correlation.R")
      plot.pairwise.correlation.heatmap(ConcordanceMatrix,
                                        cellwidth = 40,
                                        cellheight = 40,
                                        fig.width = 6,
                                        fig.height = 6,
                                        display_numbers = TRUE,
                                        outdir = workdir,
                                        fig.prefix = paste0(prefix,"_",cluster_num,"cluster_concordance")) 
      
    }
    
    # nmi_values = suppressWarnings(ConcordanceMatrix[upper.tri(ConcordanceMatrix)])
    gg = c()
    nmi_values = c()
    for(i in 1:nrow(ConcordanceMatrix)){
      for(j in 1:ncol(ConcordanceMatrix)) {
        if(j > i){
          gg = c(gg, paste0(rownames(ConcordanceMatrix)[i],"_graph and ",colnames(ConcordanceMatrix)[j],"_graph"))
          nmi_values = c(nmi_values, ConcordanceMatrix[i,j])
          # print(gg)
        }
      }
    }
    
    nmi_df = data.frame(Pair=gg,NMI=nmi_values)
    nmi_df$ClusterNum = cluster_num
    nmi_list[[cluster_num]] = nmi_df
  }
  
  nmi_all = do.call(rbind, nmi_list)
  # nmi_all$Pair = gsub("Fused_all","Fusion",nmi_all$Pair)
  write.table(nmi_all, file=file.path(workdir, "nmi_changes.tsv"), row.names = F, sep = "\t", quote = F)
  
  return(nmi_all)
}




#' Plot the NMI score changes with different numbers of clusters
#'
#' @param nmi_df 
#' @param fig.width 
#' @param fig.height 
#' @param workdir 
#' @param prefix 
#'
#' @return
#' @export
#'
#' @examples
plot_nmi_changes = function(nmi_df, 
                            fig.width = 6, 
                            fig.height = 8, 
                            workdir="./", 
                            prefix="snf"){

  # color_list = c('#636EFA','#EF553B','#00CC96',"#801000")
    
  suppressWarnings(ggplot(nmi_df, aes(x=ClusterNum, y=NMI, color=Pair)) + geom_line(size=1) + 
    ylab("Normalized Mutual Information (NMI)") +
    xlab("Number of clusters") +
    scale_x_discrete(limit = min(nmi_df$ClusterNum):max(nmi_df$ClusterNum), labels = min(nmi_df$ClusterNum):max(nmi_df$ClusterNum)) +
    ggtitle("The NMI score changes with different numbers of clusters") + 
    # scale_color_manual(values = color_list) +
    ggsci::scale_color_d3(palette = "category20") +
    ggsci::scale_fill_d3(palette = "category20") +
    theme_bw() + 
    theme(legend.position = "right",
          # legend.box="vertical", 
          # legend.margin=margin(),
          axis.title = element_text(face = "bold", size = 14),
          axis.text = element_text(size = 12, colour = "black"),
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 8),
          # panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
          # panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey"),
          # panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey")
    )) #+ guides(fill=guide_legend(nrow=2,byrow=TRUE),color=guide_legend(nrow=2,byrow=TRUE))
  
  ggsave(filename = file.path(workdir, paste0(prefix, "_nmi_changes.png")), width = fig.width, height = fig.height, dpi = 300, units = "in")
  
}




#' Plot sankeyplot for cluster concordance when incrementing the cluster numbers
#'
#' @param group_list 
#' @param workdir 
#' @param prefix 
#'
#' @return
#' @export
#'
#' @examples
plot_cluster_sankey = function(group_list, 
                               return_data = FALSE,
                               workdir="./", 
                               prefix="cluster_concordance"){
  suppressMessages(library(networkD3))
  suppressMessages(library(tidyverse))
  suppressMessages(library(plotly))

  # group_list
  cls_df = as.data.frame(do.call(cbind, group_list))
  cls_df2 = do.call(cbind,lapply(1:ncol(cls_df), function(x){
    paste(paste0("SNF",x+1),cls_df[,x],sep="_C")
  })) %>% as.data.frame()
  colnames(cls_df2) = paste("SNF",gsub("_clusters","",colnames(cls_df)),sep="_")
  
  
  # rownames(cls_df2) = smp_lst
  cls_df2$amount = rep(1,nrow(cls_df2))
  
  cls_nodes <- cls_df2 %>% 
    pivot_longer(-amount, values_to = "name_node") %>% 
    distinct(name_node) %>%
    mutate(idx = (1:n()) - 1)
  # print(cls_df2)
  cls_links = bind_rows(
    cls_df2 %>% select(source = SNF_2, target = SNF_3, amount),
    cls_df2 %>% select(source = SNF_3, target = SNF_4, amount),
    cls_df2 %>% select(source = SNF_4, target = SNF_5, amount),
    cls_df2 %>% select(source = SNF_5, target = SNF_6, amount),
    cls_df2 %>% select(source = SNF_6, target = SNF_7, amount),
    cls_df2 %>% select(source = SNF_7, target = SNF_8, amount),
    cls_df2 %>% select(source = SNF_8, target = SNF_9, amount),
    cls_df2 %>% select(source = SNF_9, target = SNF_10, amount)) %>%
    group_by(source, target) %>%
    summarise(value = sum(amount), .groups = "drop") 
  
  # cls_links <- rbind(cls_links,
  #                    data.frame(
  #                      "source" = unique(cls_df2$SNF_10),
  #                      "target" = unique(cls_df2$SNF_10),
  #                      "value" = 10**-10
  #                    ))
  # 
  # print(cls_links)
  cls_links = cls_links %>%
    mutate(across(c(source, target), ~ cls_nodes$idx[match(.x, cls_nodes$name_node)]))
  
  
  # Plot
  sankey_plot = plot_ly(
    type = "sankey",
    orientation = "h",
    node = list(label = gsub("SNF","",cls_nodes$name_node), pad = 15, thickness = 15),
    # node = list(label = cls_nodes$name_node, pad = 15, thickness = 15),
    link = as.list(cls_links))
  
  networkD3::saveNetwork(sankey_plot, file = file.path(workdir,"sankeyplot.html")) #save to html
  # webshot::webshot(file.path(workdir,"cluster_concordance_changes.sankeyplot.pdf"), vwidth = 800, vheight = 500) #save to pdf
  webshot::webshot(file.path(workdir,"sankeyplot.html"), 
                   file = file.path(workdir,paste0(prefix,"_sankey.png")), 
                   vwidth = 800, vheight = 400, zoom = 4) #save to pdf
  unlink(file.path(workdir,"sankeyplot.html"))
  unlink(file.path(workdir,"sankeyplot_files"),recursive = TRUE)
  
  if(return_data){
    return(list("links"=cls_links, "nodes"=cls_nodes))
  }
  
  
}




.discretisation = function(eigenVectors) {
  
  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors,1,normalize))
  
  n = nrow(eigenVectors)
  k = ncol(eigenVectors)
  
  R = matrix(0,k,k)
  R[,1] = t(eigenVectors[round(n/2),])
  
  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }
  
  c = matrix(0,n,1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
    i = mini(c)
    R[,j] = t(eigenVectors[i,])
  }
  
  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)
    
    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]
    
    NcutValue = 2 * (n-sum(S))
    if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps) 
      break
    
    lastObjectiveValue = NcutValue
    R = V %*% t(U)
    
  }
  
  return(list(discrete=eigenDiscrete,continuous =eigenVectors))
}  



.discretisationEigenVectorData <- function(eigenVector) {
  
  Y = matrix(0,nrow(eigenVector),ncol(eigenVector))
  maxi <- function(x) {
    i = which(x == max(x))
    return(i[1])
  }
  j = apply(eigenVector,1,maxi)
  Y[cbind(1:nrow(eigenVector),j)] = 1
  
  return(Y)
  
}



estimate_number_of_clusters_given_graph = function(W, ncluster=2:10) {
  # Estimates the best number of clusters from a vector of choices, using 
  #   the eigen-gap & rotation cost heuristics.
  #
  # Args:
  #   W: Affinity matrix (usually result from SNF)
  #   ncluster: A vector of integers specifying which cluster numbers to check
  #
  # Returns:
  #   A vector of the top two suggested number of clusters using
  #       the eigen-gap and rotation cost heuristics. 
  #
  
  #Why is this performed here?
  W <- (W + t(W))/2
  diag(W) <- 0

  # compute unnormalized Laplacian
  degs <- rowSums(W)
  degs[degs == 0] <- .Machine$double.eps    
  D <- diag(degs)    
  L <- D - W
  Di <- diag(1 / sqrt(degs))
  L <- Di %*% L %*% Di
  #print(dim(L))
  
  # compute the eigenvectors corresponding to the k smallest
  eigs <- eigen(L)
  eigs_order <- sort(eigs$values, index.return=T)$ix
  eigs$values <- eigs$values[eigs_order]
  eigs$vectors <- eigs$vectors[, eigs_order]
  eigengap <- abs(diff(eigs$values))
  #    eigengap <- eigengap * (1 - eigs$values[1:length(eigs$values) - 1]
  #        ) / (1 - eigs$values[2:length(eigs$values)])
  
  quality <- list()
  for (c_index in 1:length(ncluster)) {
    ck <- ncluster[c_index]
    UU <- eigs$vectors[, 1:ck]
    EigenvectorsDiscrete <- .discretisation(UU)[[1]]
    EigenVectors <- EigenvectorsDiscrete^2
    
    #MATLAB: sort(EigenVectors,2, 'descend');
    temp1 <- EigenVectors[do.call(order, lapply(1:ncol(EigenVectors),
                                                function(i) EigenVectors[, i])), ]
    temp1 <- t(apply(temp1, 1, sort, TRUE))  
    
    quality[[c_index]] <- (1 - eigs$values[ck + 1]) / 
      (1 - eigs$values[ck]) * 
      sum( sum( diag(1 / (temp1[, 1] + .Machine$double.eps) ) %*%
                  temp1[, 1:max(2, ck-1)] ))
  }
  #Eigen-gap best two clusters
  m1 <- eigengap[ncluster]
  t1 <- sort(eigengap[ncluster], decreasing=TRUE, index.return=T)$ix
  K1 <- ncluster[t1[1]]
  K12 <- ncluster[t1[2]]
  
  #Rotation cost best two clusters
  m2 <- quality
  t2 <- sort(unlist(quality), index.return=TRUE)$ix
  K2 <- ncluster[t2[1]]
  K22 <- ncluster[t2[2]]    
  
  output <- list("Eigen-gap best"=K1, "Eigen-gap 2nd best"=K12,
                 "Rotation cost best"=K2, "Rotation cost 2nd best"=K22)
  output_2 <- list(m1, m2, K1, K12, K2, K22)
  return (output_2)
}





#' Using label propagation to predict the subtypes/labels of new patients/data points - predicting the new labels with label propagation.
#'
#' @param dat Data list, new patients must be added into the data list along with existed patients, and have the same number of omic view  
#' @param K Number of neighbors.
#' @param alpha Hyperparameter used in constructing similarity network.
#' @param method A indicator of which method to use to predict the label. method = 0 means to use local and global consistency; method = 1 means to use label propagation.
#' @param t Number of iterations.
#'
#' @return The prediction of which group the new patients/data points belongs to.
#' @export
#'
#' @examples
do_snf_propagation2predict = function(dataL=list(),
                                      label=NULL,
                                      K=20,
                                      alpha=0.5,
                                      t=20,
                                      method=TRUE){
  
  # Predicting the new labels with label propagation
  # Note: The training data must have the same number of view and columns as the test data?
 
  # Create the training and test data
  n = floor(0.8*length(label)) # number of training cases
  trainSample = sample.int(length(label), n)
  train = lapply(dataL, function(x) x[trainSample, ]) # Use the first 150 samples for training
  test = lapply(dataL, function(x) x[-trainSample, ]) # Test the rest of the data set
  groups = label[trainSample]
  
  # Apply the prediction function to the data
  newLabel = SNFtool::groupPredict(train, test, groups, K, alpha, t, method)
  
  # Compare the prediction accuracy
  accuracy = sum(label[-trainSample] == newLabel[-c(1:n)])/(length(label) - n)
  
  return(list(newLabel=newLabel,
              accuracy=accuracy))
}



do_exectueSNF = function(dat,
                         clusterNum=3,
                         idx=4,
                         workdir="./",
                         prefix="exectueSNF"){
  
  print(paste0("*** working on multi-omics data using cluster number equal to ",clusterNum))
  result = CancerSubtypes::ExecuteSNF(dat,
                                      clusterNum = clusterNum,
                                      K=20,
                                      alpha=0.5,
                                      t=20,
                                      plot = TRUE)
  
  # save(result, file = file.path(workdir,paste0("integrated_multiomics_",clusterNum,"cluster_result.RData")))
  # print(paste0("*** save ",level," result for clusterNum=",clusterNum))
  # 
  # load(file.path(workdir,paste0("integrated_multiomics_",clusterNum,"cluster_result.RData")))
  # idx = 4
  # sil = CancerSubtypes::silhouette_SimilarityMatrix(result$originalResult[[idx]]$consensusClass,
  #                                                   result$originalResult[[idx]]$consensusMatrix)
  # sil_plot = plot_silhouette(sil)
  # 
  # # sil_plot
  # 
  # ggsave(sil_plot,
  #        file = file.path(workdir, paste0("integrated_multiomics_sil_plot", idx, ".png")),
  #        width = 7,
  #        height = 7,
  #        units = "in",
  #        dpi = 300)
  # 
  # 
  # setwd(curr_dir)
}




do_kernel_nmf_clustering = function(matrix_x, 
                                    matrix_y=NULL, 
                                    sigma=0.05, 
                                    C=4, 
                                    plot=FALSE, 
                                    workdir="./", 
                                    prefix="knmf") {
  
  suppressMessages(library(kernlab))
  suppressMessages(library(NMF))
  
  ## create a RBF kernel function with sigma hyper-parameter 0.05 
  rbf = rbfdot(sigma = sigma)
  
  ## Compute the kernel matrix using a radial basis function (RBF) kernel
  if(is.null(matrix_y)){
    kernel_matrix = kernelMatrix(rbf, x = matrix_x)
  }else{
    kernel_matrix = kernelMatrix(rbf, x = matrix_x, y = matrix_y)
  }
  
  # Perform KNMF using the kernel matrix
  knmf_result = nmf(kernel_matrix, rank = C)
  
  # Get the factorization matrices
  W_basis = basis(knmf_result)
  H_coef = coef(knmf_result)
  
  if(plot){
    pdf(file=file.path(workdir,paste0(prefix,"_knmf_",C,"_cluster.pdf")))
    # res = coefmap(knmf_result)
    res = consensusmap(knmf_result,color = "-RdYlBu")
    dev.off()
  }
  
  # # get the final subtypes information
  # group = 	# the final subtypes information
  
  return(knmf_result)
}














                             









