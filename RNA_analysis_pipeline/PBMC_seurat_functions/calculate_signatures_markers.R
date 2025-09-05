source("/home/tferrari/NAS2/iris/1_scripts/iris_scRNAseq/2_functions/scRNA-seq/PBMC_seurat_functions/create_and_save_plots.R")

calculate_cell_signatures <- function(signatures, seurat.object, ngenes, reduction.name, path){
  signatures$cluster <- gsub(" ", "_", signatures$cluster)
  
  signatures %>% group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.25) %>%
    slice_head(n = ngenes) %>%
    ungroup() -> top.PBMC.DEGs

    #Split the table according to clusters
    l.signature.PBMCs <- split(top.PBMC.DEGs, top.PBMC.DEGs$cluster)
    lapply(l.signature.PBMCs, nrow)

    #Score each cell for cell type signatures
    RNA.norm <- as.data.frame(GetAssayData(seurat.object, assay = "RNA"))
    scaled <- apply(RNA.norm, 1, scale)
    rownames(scaled) <- colnames(RNA.norm)
    scaled <- t(scaled)

    for(i in seq_along(l.signature.PBMCs)) {
      signature.i <- l.signature.PBMCs[[i]]$gene
      subset.i <- names(l.signature.PBMCs)[i]
      print(subset.i)

      #Calculate and store zscore
      scaled.sig <- subset(scaled, rownames(scaled) %in% signature.i)

      #Check for NAs and switch their values to 0 if necessary
      scaled.sig.per.cell <- colSums(scaled.sig, na.rm = TRUE) / nrow(scaled.sig)
      scaled.sig.per.cell[is.nan(scaled.sig.per.cell)] <- 0
      # scaled.sig.per.cell[scaled.sig.per.cell <= 0] <- 0
      
      print(summary(scaled.sig.per.cell))

      seurat.object@meta.data[, subset.i] <- as.numeric(scaled.sig.per.cell)
    }
    
    if(!file.exists(path)){
      dir.create(path, recursive = TRUE)
      message(paste("Directory created: ", path))
    } else {
      message(paste("Directory already exists: ", path))
    }

    #Plot the feature per signature cluster group
    for(i in seq_along(names(l.signature.PBMCs))){
    p <- FeaturePlot(object = seurat.object,
                     features = names(l.signature.PBMCs)[i],
                     reduction = reduction.name,
                     cols = c("lightgrey", "deepskyblue", "blue"),
                     pt.size = 0.8
                     ) +
      labs(title = NULL) +
      plot_annotation(title = paste("Z-score projection for", 
                                    gsub("_", " ", names(l.signature.PBMCs[i])))) &
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    print(p)
    ggsave(filename = paste0(path, "/ftr_plot_", i, ".png"), plot = p, width = 8, height = 8)
    }
    
    return(seurat.object)
}

calculate_and_plot_marker_DEGs <- function(seurat.object, min.pct, logfc.threshold,
                                           num.mrkrs, path){
  markers.all <- FindAllMarkers(seurat.object, only.pos = FALSE,
                                min.pct = min.pct, logfc.threshold = logfc.threshold)
  
  markers.all <- markers.all %>% arrange(cluster, desc(avg_log2FC))
  head(markers.all)

  #Select top 250 DEGs
  markers.all %>%
    group_by(cluster) %>%
    top_n(n = num.mrkrs, wt = avg_log2FC) -> topn
  #Select restricted top DEGs for plotting
  #markers.all %>% group_by(cluster) %>%
  #  top_n(n = 20, wt = avg_log2FC) -> top20
  #markers.all %>% group_by(cluster) %>%
  #  top_n(n = 10, wt = avg_log2FC) -> top10
  markers.all %>% group_by(cluster) %>%
    top_n(n = 7, wt = avg_log2FC) -> top7

  markers.all %>%
    group_by(cluster) %>%
    slice_max(n=num.mrkrs, order = -avg_log2FC) -> bottomn

  #Plot heatmap
  n_classes = length(unique(Idents(seurat.object)))
  p <- DoHeatmap(seurat.object, features = top7$gene,
              size = 3, angle = 90) + #,
      #         group.colors = c("deepskyblue","red2","green","green3","red","deepskyblue3", "grey60", "grey80")) + 
    scale_fill_gradient2( low = rev(c('#436EEE','#3A5FCD','#27408B')),
                          mid = "white", high = rev(c('#FFD700','#EEC900','#CDAD00')),
                          midpoint = 0, guide = "colourbar", aesthetics = "fill") +
    NoLegend()
  print(p)
  ggsave(filename = path, plot = p, width = n_classes+8, 
         height = n_classes+8)

  return(list(markers.all, topn, bottomn))
}