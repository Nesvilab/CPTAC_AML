source("../../../cluster/silhouette_analysis.R")
library(umap)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)



ratio_sil_changes = data.frame()


n_componements = 3
slice_ratios = seq(0.1, 1, by=0.1)


for(ratio in slice_ratios){
  sil_changes = c()

  snf_result = readRDS(file.path("D:/Projects/CPTAC/AML/AML_Data/Figures/rna_protein_snf_clustering_update_correction_alec2",
                                 paste0("individual/RNA.Protein-",ratio,"/RNA.Protein-",ratio,"-individual.RDS")))

  W = snf_result[["dat_graph"]][["fusion_graph"]]
  custom.settings = umap.defaults
  custom.settings$n_neighbors = 20
  w_umap = umap::umap(W, n_components = n_componements, random_state = 15, config = custom.settings)

  for(C in 2:10){
    group = snf_result[["group"]][[paste0(C,"_clusters")]]
    res_df = as.data.frame(w_umap[["layout"]])
    res_df$ID = colnames(W)
    res_df$Cluster = group

    if(!dir.exists(file.path(workdir, paste0("RNA.Protein-",ratio)))) dir.create(file.path(workdir, paste0("RNA.Protein-",ratio)), recursive = T)
    write.table(res_df,
                file = file.path(workdir, paste0("RNA.Protein-",ratio,"/umap",n_componements,"_",C,"clusters.tsv")),
                sep="\t",
                quote = F,
                row.names = F)

    umap_result = calculate_umap_dist(
      filename = file.path(workdir, paste0("RNA.Protein-",ratio,"/umap",n_componements,"_",C,"clusters.tsv"))
      )
    saveRDS(umap_result,
            file = file.path(workdir,
                             paste0("RNA.Protein-",ratio,"/umap",n_componements,"_",C,"clusters.RDS")))

    umap_dist = umap_result[["dist"]]
    umap_group = umap_result[["group"]]
    sil_result = plot_silhouette(silhouette(umap_group, umap_dist))
    sil_plot = sil_result[["plot"]]
    sil_df = sil_result[["sil_df"]]
    write.table(sil_df,
                file = file.path(workdir,paste0("sil_umap",n_componements,"_",C,"clusters.tsv")),
                sep = "\t",
                quote = F,
                row.names = F)

    # sil_plot
    ggsave(sil_plot,
           filename=file.path(workdir,
                              paste0("sil_umap",n_componements,"_",ratio,"_",C,"clusters.png")),
           width = 8,
           height = 8,
           dpi = 300,
           units = "in")

    sil_changes = c(sil_changes, round(mean(sil_df$sil_width), 2))

  }

  if(nrow(ratio_sil_changes)==0){
    ratio_sil_changes = data.frame(matrix(sil_changes, nrow = 1))
  }else{
    ratio_sil_changes = rbind(ratio_sil_changes,data.frame(matrix(sil_changes, nrow = 1)))
  }
}

colnames(ratio_sil_changes) = 2:10
rownames(ratio_sil_changes) <- paste0(as.numeric(slice_ratios)*100,"%", " input data")
colnames(ratio_sil_changes) <- paste0(as.numeric(colnames(ratio_sil_changes), " clusters"))


rownames(ratio_sil_changes) <- gsub(" input data","",rownames(ratio_sil_changes))
ratio_sil_changes <- ratio_sil_changes[,-1]

ht_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(50)
ht <- Heatmap(
  as.matrix(ratio_sil_changes),
  name = "Value",
  cluster_rows = F,
  cluster_columns = F,
  col = ht_colors,
  show_column_names = FALSE,
  row_names_side = "left",
  bottom_annotation = HeatmapAnnotation(
    text = anno_text(
      colnames(ratio_sil_changes),
      rot = 0,
      offset = unit(0.5, "npc"),
      just = "center"
    ),
    annotation_height = max_text_width(colnames(ratio_sil_changes))
  ), 
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(ratio_sil_changes[i, j], x, y, gp = gpar(fontsize = 10))
  },
  rect_gp = gpar(col = "white", lwd = 0.5))




pdf(file.path(workdir, paste0("sil_umap",n_componements,"_changes.pdf")), 
    height = 3.8, width = 4.2)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.95, height=0.95, name="vp", just=c("right","top"))), action="prepend")
draw(ht, show_heatmap_legend = TRUE, column_title="Mean silhouette width changes")

setHook("grid.newpage", NULL, "replace")
grid.text("Examined cluster number", y=-0.01, gp=gpar(fontsize=12))
grid.text("The percentage of input data", x=-0.01, rot=90, gp=gpar(fontsize=12))

dev.off()