suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(ggsci))
suppressMessages(library(survminer))



plot_silhouette = function(sil, color = "red") {
  
  temp_data = data.frame(cluster = sil[, 1],
                         neighbor = sil[, 2],
                         sil_width = sil[, 3],
                         stringsAsFactors = FALSE)
  
  temp_data = temp_data %>% 
    dplyr::mutate(cluster = as.character(cluster)) %>%
    dplyr::arrange(desc(cluster), sil_width) %>%
    dplyr::mutate(index = 1:nrow(temp_data))
  
  plot = temp_data %>% 
    ggplot() +
    geom_bar(aes(y = sil_width,x = index,fill = cluster,color = cluster),
             stat = "identity",
             show.legend = FALSE
    ) +
    geom_hline(yintercept = 0) +
    ggsci::scale_color_npg() +
    ggsci::scale_fill_npg() +
    scale_y_continuous(expand = expansion(mult = c(0, .2))) +
    theme_classic() +
    labs(y = paste("Silhouette width","\nAverage silhouettle width:",round(mean(temp_data$sil_width), 2)), x = "") +
    theme(
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 13),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank()) +
    coord_flip()
  
  plot = plot +
    ggplot2::annotate(geom = "text",
                      y = 0,
                      x = max(temp_data$index),
                      label = paste("n =", nrow(temp_data)),
                      color = "black",
                      hjust = -0.5,
                      vjust = -1,
                      size = 4)
  
  
  cluster_num = as.numeric(max(temp_data$cluster))
  title = paste(cluster_num, "clusters", "\n", "clusterIndex: size | Sil_width")
  
  plot = plot +
    ggplot2::annotate(geom = "text",
                      y = max(temp_data$sil_width),
                      x = max(temp_data$index),
                      label = title,
                      color = "black",
                      hjust = 0,
                      vjust = 0,
                      size = 4)
  
  class = temp_data$cluster %>%
    unique() %>%
    sort() %>%
    rev()
  
  cluster_num = temp_data %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(cluster)) %>%
    dplyr::pull(n)
  
  cluster_num = sapply(1:length(cluster_num), function(x) {
    if (x == 1) {
      cluster_num[x] / 2
    } else{
      tail(cumsum(cluster_num[1:(x - 1)]), 1)  +  cluster_num[x] / 2
    }
  })
  
  for (i in 1:length(class)) {
    label = paste(class[i],":",sum(temp_data$cluster == class[i]),"|",
                  round(mean(temp_data$sil_width[temp_data$cluster == class[i]]), 2))
    plot = plot +
      ggplot2::annotate(geom = "text",
                        y = max(temp_data$sil_width),
                        x = cluster_num[i],
                        label = label,
                        color = "black",
                        hjust = 0,
                        vjust = 0,
                        size = 4)
    
  }
  
  return(list(plot=plot,
              sil_df=temp_data)) 
}



plot_silhouette_changes = function(sil_df,
                                   fig.width = 6, 
                                   fig.height = 4, 
                                   workdir="./", 
                                   prefix="sil"){
  
  suppressWarnings(ggplot(sil_df, aes(x=clusterNum, y=sil_mean_width)) + 
                     geom_line(size=1, color="dodgerblue") + 
                     ylab("Mean silhouette score") +
                     xlab("Number of clusters") +
                     ylim(0,1) + 
                     scale_x_discrete(limit = min(sil_df$clusterNum):max(sil_df$clusterNum), labels = min(sil_df$clusterNum):max(sil_df$clusterNum)) +
                     ggtitle("The silhouette score changes with different numbers of clusters") + 
                     # scale_color_manual(values = color_list) +
                     ggsci::scale_color_d3(palette = "category20") +
                     ggsci::scale_fill_d3(palette = "category20") +
                     theme(axis.title = element_text(face = "bold", size = 14),
                           axis.text = element_text(size = 12, colour = "black"),
                           legend.title = element_text(face = "bold", size = 14),
                           legend.text = element_text(size = 12),
                           panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
                           # panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey"),
                           # panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey")
                     ))
  
  ggsave(filename = file.path(workdir, paste0(prefix, "_sil_changes.png")), width = fig.width, height = fig.height, dpi = 300, units = "in")
  
}


plot_umap = function(W, 
                     C=10,
                     group=NULL, 
                     n_components=2, 
                     n_neighbors = 20,
                     use_centroids=TRUE,
                     add_ellipse=FALSE, 
                     add_segments=FALSE,
                     workdir="./",
                     prefix="umap"){
  
  suppressMessages(library(umap))
  suppressMessages(library(plotly))
  
  # w_umap = umap::umap(W, n_components = 2, random_state = 15)
  # layout = w_umap[["layout"]]
  # layout = data.frame(layout)
  # gg = cbind(layout, Cluster=paste0("C",group))
  # 
  # plt = plot_ly(gg, x = ~X1, y = ~X2, split = ~Cluster, colors = ggsci::pal_jco()(10), type = 'scatter', mode = 'markers', size = 5)%>%
  #   layout(
  #     plot_bgcolor = "#e5ecf6",
  #     # legend=list(title=list(text='Cluster')),
  #     xaxis = list(
  #       title = "UMAP-1"),
  #     yaxis = list(
  #       title = "UMAP-2"))
  # plt
  if(n_components==2){
    custom.settings = umap.defaults
    custom.settings$n_neighbors = n_neighbors
    w_umap = umap::umap(W, n_components = 2, random_state = 15, config = custom.settings) 
    res_df = as.data.frame(w_umap[["layout"]])
    colnames(res_df) = c("UMAP1","UMAP2")
    res_df$ID =rownames(W)
    res_df$group = group
    # res_df$group = paste0("C",group)
    
    if(C<5){
      col_list = c('#636EFA','#EF553B','#00CC96',"#801000") # c("#2273C3","#EFC144","#868686")
    }else{
      col_list = ggsci::pal_jco()(20)
    }
    
    gg = data.frame(Cluster=factor(res_df$group), x=res_df$UMAP1, y=res_df$UMAP2, ID=res_df$ID)
    # rownames(gg) = res_df$ID
    centroids =  aggregate(cbind(x,y)~Cluster, data=gg, mean)
    gg = merge(gg,centroids,by="Cluster",suffixes=c("",".centroid"))
    write.table(gg, file.path(workdir,paste0(prefix,"_umap.tsv")), row.names = F, sep="\t",quote = F)
    
    plt = ggplot(gg)+
      geom_point(aes(x = x, y = y, color = Cluster, fill = Cluster), size=3)+
      # ggthemes::theme_base() +
      theme_bw() +
      theme(legend.title=element_text(size=15),
            legend.text=element_text(size=12),
            legend.key.size=unit(0.7,"line"),
            axis.title = element_text(size=16, face = "bold",color="black"),
            axis.text = element_text(color="black",size=16, face = "bold"),
            strip.text.x = element_text(size = 9,color="black"),
            plot.title = element_text(hjust = 0.5, size = 20, face = "bold", vjust=2),
            legend.background=element_rect(fill = "white"))+ 
      xlab("UMAP1") +
      ylab("UMAP2") +
      scale_color_manual(values=col_list) + 
      scale_fill_manual(values=col_list)
    
    if(use_centroids){
      plt = plt + geom_point(data=centroids, aes(x=x, y=y, color=Cluster, fill = Cluster), size = 5, alpha = 0.9, shape = 19)
    }
    if(add_ellipse){
      plt = plt + stat_ellipse(geom = "polygon", data = gg, aes(x = x, y = y, color = Cluster, fill = Cluster), alpha = 0.3) 
    }
    if(add_segments){
      plt = plt + geom_segment(data = gg, aes(x=x.centroid, y=y.centroid, xend=x, yend=y, color=Cluster)) +
        annotate("text",x=gg$x.centroid-0.2,y=gg$y.centroid+0.2,label=gg$Cluster)
    }
    
    plt = ggplotly(plt)
    plt$x$data = lapply(plt$x$data, function(x){
      if(x$mode=="lines"){
        x$showlegend = FALSE
      }
      return(x)
    })
    # save_image(plt, file = file.path(workdir,"umap.png"), width = 8 * 96, height = 8 * 96)
    # plotly::export(p = plt, #the graph to export
    #                file = file.path(workdir,"umap.png")) #the name and type of file (can be .png, .jpeg, etc.)
    
    # ggsave(filename = file.path(workdir,"umap.png"), width = 8, height = 8, dpi=300, units = "in")
    
    return(plt)
    
  }else if(n_components == 3){
    w_umap = umap::umap(W, n_components = 3, random_state = 15) 
    res_df = as.data.frame(w_umap[["layout"]])
    colnames(res_df) = c("UMAP1","UMAP2","UMAP3")
    res_df$group = paste0("C",group)
    
    gg = data.frame(Cluster=factor(res_df$group), res_df[,c("UMAP1","UMAP2","UMAP3")])
    
    if(C<5){
      col_list <- c('#636EFA','#EF553B','#00CC96',"#801000") # c("#2273C3","#EFC144","#868686")
    }else{
      col_list = ggsci::pal_jco()(10)
    }
    
    plt2 = plot_ly(gg, x = ~UMAP1, y = ~UMAP2, z = ~UMAP3, split = ~Cluster, colors =  col_list[1:length(unique(group))]) 
    plt2 = plt2 %>% add_markers() 
    plt2 = plt2 %>% layout(scene = list(xaxis = list(title = 'UMAP1'),
                                        yaxis = list(title = 'UMAP2'),
                                        zaxis = list(title = 'UMAP3'))) 
    
    return(plt2)
  }
  
}


calculate_umap_dist = function(filename=NULL){
  euclidean_dist = function(x, y) sqrt(sum((x - y)^2))
  
  tmp = read.delim(filename, sep="\t")
  
  # init umap_dist matrix
  umap_dist = data.frame(matrix(0,nrow = nrow(tmp), ncol = nrow(tmp)))
  rownames(umap_dist) = sort(tmp$ID)
  colnames(umap_dist) = sort(tmp$ID)
  
  for(r in rownames(umap_dist)){
    # smp1 = tmp[tmp$ID==r,c("x","y")]
    smp1 = tmp[tmp$ID==r, grepl("^V",colnames(tmp))]
    for(c in colnames(umap_dist)){
      smp2 = tmp[tmp$ID==c, grepl("^V",colnames(tmp))]
      umap_dist[r,c] = euclidean_dist(smp1,smp2)
    }
  }
  # print(head(umap_dist)[,1:5])
  group = tmp[match(rownames(umap_dist),tmp$ID),]$Cluster
  return(list(dist=umap_dist,group=group))
}








