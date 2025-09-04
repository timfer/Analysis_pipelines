# generate_color_palette <- function(n_groups) {
#   # Define the base colors for the first 8 discrete groups
#   base_colors <- c("deepskyblue", "red2", "green", "green3", "red", 
#                    "deepskyblue3", "grey60", "grey80")
#   
#   # If groups <= 8, return the base colors
#   if (n_groups <= length(base_colors)) {
#     return(base_colors[1:n_groups])
#   } else {
#     # If groups > 8, generate additional colors with gradient transition
#     # Define gradient colors based on your `scale_fill_gradient2` scheme but without white
#     gradient_colors <- colorRampPalette(rev(c('#d1e5f0', '#67a9cf', '#2166ac', 
#                                               '#b2182b', '#ef8a62', '#fddbc7')))(n_groups - length(base_colors))
#     
#     # Combine base colors and gradient colors
#     return(c(base_colors, gradient_colors))
#   }
# }

generate_QC_plots <- function(seurat.object, res, reduction.name, target_expID, path) {
  # Ensure the output directory exists
  if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
  }
  # n_clusters <- length(seurat.object@meta.data$seurat_clusters)
  # cluster_colors <- generate_color_palette(n_clusters)
  # 
  p <- VlnPlot(seurat.object, features = c("nFeature_RNA"),
               group.by = "expID", pt.size = 0) + 
    labs(title = NULL) +
    xlab("Experimental ID (batch)") + ylab("nFeature_RNA (counts)") + 
    plot_annotation(title = "Feature RNA counts per experimental batch") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    NoLegend() +
    geom_boxplot(width = 0.1, fill = NA, color = "black", position = position_dodge(width = 0.9))
  print(p)
  ggsave(filename = paste0(path, "/nftr_RNA_vln_expID.png"), plot = p,
         width = 8, height = 8)
  
  p <- VlnPlot(seurat.object, features = c("nCount_RNA"),
               group.by = "expID", pt.size = 0) + 
    labs(title = NULL) +
    xlab("Experimental ID (batch)") + ylab("nCount_RNA (counts)") + 
    plot_annotation(title = "RNA counts per experimental batch") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    NoLegend() +
    geom_boxplot(width = 0.1, fill = NA, color = "black", position = position_dodge(width = 0.9)) +
    scale_y_log10(breaks = 10^(1:6), labels = comma(10^(1:6)))
  print(p)
  ggsave(filename = paste0(path, "/ncount_RNA_vln_expID.png"), plot = p,
         width = 8, height = 8)
  
  p <- VlnPlot(seurat.object, features = c("percent.mt"),
               group.by = "expID", pt.size = 0) + 
    labs(title = NULL) +
    xlab("Experimental ID (batch)") + ylab("Mitochondrial gene expression (per cent (%))") + 
    plot_annotation(title = "Fraction of mitochondrial genes expressed per experimental batch") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    NoLegend() +
    geom_boxplot(width = 0.1, fill = NA, color = "black", position = position_dodge(width = 0.9))
  print(p)
  ggsave(filename = paste0(path, "/pct_mt_vln_expID.png"), plot = p,
         width = 8, height = 8)
  
  p <- VlnPlot(seurat.object, features = c("nFeature_RNA"),
               group.by = "seurat_clusters", pt.size = 0) + 
    labs(title = NULL) +
    xlab(paste0("UMAP cluster (resolution ", res,")")) + ylab("nFeature_RNA (counts)") + 
    plot_annotation(title = "Feature RNA counts per integrated low-dimensional cluster") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    NoLegend() +
    geom_boxplot(width = 0.1, fill = NA, color = "black", position = position_dodge(width = 0.9))
  print(p)
  ggsave(filename = paste0(path, "/nftr_RNA_vln_cluster.png"), plot = p,
         width = 8, height = 8)
  
  p <- VlnPlot(seurat.object, features = c("nCount_RNA"),
               group.by = "seurat_clusters", pt.size = 0) + 
    labs(title = NULL) +
    xlab(paste0("UMAP cluster (resolution ", res,")")) + ylab("nCount_RNA (counts)") + 
    plot_annotation(title = "RNA counts per integrated low-dimensional cluster") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    NoLegend() +
    geom_boxplot(width = 0.1, fill = NA, color = "black", position = position_dodge(width = 0.9))
  print(p)
  ggsave(filename = paste0(path, "/ncount_RNA_vln_cluster.png"), plot = p,
         width = 8, height = 8)
  
  p <- VlnPlot(seurat.object, features = c("percent.mt"),
               group.by = "seurat_clusters", pt.size = 0) + 
    labs(title = NULL) +
    xlab(paste0("UMAP cluster (resolution ", res,")")) + ylab("Mitochondrial gene expression (per cent (%))") + 
    plot_annotation(title = "Feature RNA counts per integrated low-dimensional cluster") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    NoLegend() +
    geom_boxplot(width = 0.1, fill = NA, color = "black", position = position_dodge(width = 0.9))
  print(p)
  ggsave(filename = paste0(path, "/pct_mt_vln_cluster.png"), plot = p,
         width = 8, height = 8)
  
  p1 <- FeatureScatter(seurat.object, feature1 = "n_mapped_reads", feature2 = "n_UMI",
                       group.by = "expID") # +
    # scale_x_continuous(labels = function(x) x / 1e4) + xlab(expression("n_mapped_reads (×10"^4*")")) +
    # scale_y_continuous(labels = function(x) x / 1e2) + ylab(expression("n_UMI (×10"^2*")"))
  p2 <- FeatureScatter(seurat.object, feature1 = "n_mapped_reads", feature2 = "nFeature_RNA",
                       group.by = "expID") # +
    # scale_x_continuous(labels = function(x) x / 1e4) + xlab(expression("n_mapped_reads (×10"^4*")"))
  p <- p1 + p2 + plot_layout(guides = 'collect')
  print(p)
  ggsave(filename = paste0(path, "/reads_ftr_scatter.png"), plot = p,
         width = 8, height = 8)
  
  p1 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA",
                      feature2 = "percent.mt", group.by = "expID")
  p2 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA",
                      feature2 = "nFeature_RNA", group.by = "expID")
  p <- p1 + p2 + plot_layout(guides = 'collect')
  print(p)
  ggsave(filename = paste0(path, "/var_ftr_scatter.png"), plot = p,
         width = 8, height = 8)
  
  #Plot regression variables on feature scatter plot
  ftrs = c("nCount_RNA","nFeature_RNA", "percent.mt", "percent.Hsp")
  p <- FeaturePlot(seurat.object, features = ftrs, reduction = reduction.name, 
                  pt.size = 0.3)
  print(p)
  ggsave(filename = paste0(path, "/ftr_plot_regression_vars.png"), plot = p,
         width = 16, height = 16)
  
  #Plot the UMAPs split by experimental ID
  p <- DimPlot(seurat.object, reduction = reduction.name, label = FALSE,
              pt.size = 0.5, split.by = "expID", ncol = 4) + #, cols = cluster_colors) + 
    NoLegend()
  print(p)
  ggsave(
    filename = paste0(path, "/UMAP_unnamed_exp_split.png"), 
    plot = p,
    width = 8 * min(length(unique(seurat.object@meta.data$expID)), 4), 
    height = ceiling(length(unique(seurat.object@meta.data$expID)) / 4) * 8, 
    limitsize = FALSE
  )
  
  # #Plot a single target experiment vs the rest
  # if(any(target_expID != "")){
  #   plots_list <- list()
  #   for(expID_i in target_expID){
  #     target_subset <- subset(seurat.object, expID == expID_i)
  #     p_target <- DimPlot(target_subset, reduction = reduction.name,
  #                         label = FALSE, pt.size = 0.3) + #, cols = cluster_colors) +
  #       NoLegend() + ggtitle(paste("Experiment: ", expID_i))
  #     plots_list[[length(plots_list)+1]] <- p_target
  #     }
  #   others_subset <- subset(seurat.object, cells = which(!(seurat.object@meta.data$expID %in% target_expID)))
  #   p_others <- DimPlot(others_subset, reduction = reduction.name,
  #                       label = FALSE, pt.size = 0.3) + #, cols = cluster_colors) +
  #     NoLegend() + ggtitle("Other experiments, grouped")
  #   plots_list[[length(plots_list)+1]] <- p_others
  # 
  # library(patchwork)
  # p <- wrap_plots(plots_list, ncol = length(plots_list))
  # 
  # if(length(target_expID) == 1){
  #   width <- 16
  # } else {
  #   width <- (length(target_expID)+1)*8
  # }
  # ggsave(filename = paste0(path, "/UMAP_target_exp_v_rest.png"), plot = p,
  #          width = width, height = 8, limitsize = FALSE)
  # } else {
  #   message("target_expID is empty, skipping the plot.")
  # }
  #Plot UMAP with overlayed experiments to check for batch effects
  colors <- brewer.pal(length(unique(seurat.object$expID)), "Set2") #Set colors
  p <- DimPlot(seurat.object, reduction = reduction.name, pt.size = 0.5,
              group.by = "expID")
  print(p)
  ggsave(filename = paste0(path, "/UMAP_unnamed_exp_overlap.png"), plot = p,
         width = 12, height = 8)
}

generate_bar_plots <- function(seurat.object, path) {
  # Prepare data by grouping based on Seurat cluster, 'orig.ident' (experimental ID), and 'CC'
  cluster_counts <- seurat.object@meta.data %>%
    group_by(seurat_clusters, orig.ident, CC) %>%
    tally() %>% ungroup()
  
  # Absolute Bar Plot
  p <- ggplot(cluster_counts, aes(x = seurat_clusters, y = n, fill = CC)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    facet_wrap(~ orig.ident) +
    labs(title = "Absolute Number of Cells by Cluster", 
         y = "Cell Count (absolute)", 
         x = "Seurat Cluster") # +
    # theme_minimal() +
    # theme(panel.background = element_rect(fill = "white", color = NA))
  
  print(p)
  ggsave(filename = paste0(path, "/absolute_distribution_bar_plot.png"), plot = p,
         width = 12, height = 8)
  
  cluster_counts <- cluster_counts %>%
    group_by(orig.ident, CC) %>%
    mutate(relative_count = n / sum(n)) %>%
    ungroup()
  
  # Relative Bar Plot
  p <- ggplot(cluster_counts, aes(x = seurat_clusters, y = relative_count, fill = CC)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    facet_wrap(~ orig.ident) +
    labs(title = "Relative Proportion of Cells by Cluster", 
         y = "Relative Cell Count (%)", 
         x = "Seurat Cluster") # +
    # theme_minimal() +
    # theme(panel.background = element_rect(fill = "white", color = NA))
  
  print(p)
  ggsave(filename = paste0(path, "/relative_distribution_bar_plot.png"), plot = p,
         width = 12, height = 8)
}

generate_ftr_plots_umap <- function(seurat.object, reduction.name, path){
  # Ensure the output directory exists
  if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
  }
  # n_clusters <- length(seurat.object@meta.data$seurat_clusters)
  # cluster_colors <- generate_color_palette(n_clusters)

  #####
  #Feature plots & UMAP
  #####
  #Plot cell type discriminating genes on feature plot
  ftrs = c("CD14", "CD3E", "CD19", "NKG7", "CD4",
          "TNF", "IL6", "IL1B", "CD8A", "GZMA", "IL2RA", "NCAM1", 
          "CD34", "MKI67")
  p <- FeaturePlot(seurat.object, features = ftrs, reduction = reduction.name,
                  pt.size = 0.3)
  print(p)
  ggsave(filename = paste0(path, "/ftr_plots_genes.png"), plot = p, width = 32, height = 24)
  
  #####
  #UMAPs
  #####
  # Plot the combined UMAP
  # n_clusters <- length(unique(seurat.object@meta.data$seurat_clusters))
  # umap_colors <- generate_color_palette(n_clusters)
  
  p <- DimPlot(seurat.object, label = TRUE) #, cols = umap_colors)
  print(p)
  ggsave(filename = paste0(path, "/UMAP_unnamed.png"), plot = p,
         height = 8, width = 10)
}

plot_specific_ftr <- function(seurat.object, path){
    ftr_plot_path = path
    if (!dir.exists(ftr_plot_path)) {
            dir.create(ftr_plot_path, recursive = TRUE)
    }

    #Monocyte feature plot: Done
    p <- FeaturePlot(seurat.object, features = c("HLA.DRA", "CD14", "CSF3R", "CCR2", #Classical monocytes
                                                "LYZ", "S100A8", "S100A9", "ITGAM",#Classical monocytes
                                                "FCGR3A", "CX3CR1", "CD86", "ITGB2", "ITGAL")) #Non-classical monocytes
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/monocyte_ftr_plot.png"), 
            plot = p)

    #Naive CD4 T-cell feature plot: Done
    p <- FeaturePlot(seurat.object, features = c("CD4", "LEF1", "TCF7", "CCR7", "SELL", 
                                                "LTB", "IL7R", "CD28", "FOXP1",
                                                "CD3D", "CD69", "BCL2")) 
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/naive_cd4_ftr_plot.png"), 
            plot = p)

    #Central memory CD4+ T cells
    p <- FeaturePlot(seurat.object, features = c("CD4", "IL2RA", "CD44", "SELL"))
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/central_mem_cd4_ftr_plot.png"), 
            plot = p)

    #Th2 effector memory CD4 feature plot: Done
    p <- FeaturePlot(seurat.object, features = c("GATA3", "STAT6", "STAT5A", "CCR4",
                                                "RELA", "CD40LG", "ICOS", "TNFRSF4"))
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/memory_cd4_Th2_ftr_plot.png"), 
            plot = p)

    #T-reg CD4 feature plot: Done
    p <- FeaturePlot(seurat.object, features = c("FOXP3", "IL2RA", "CTLA4", 
                                                "IKZF2", "TGFB1", "ITGB1",
                                                "LAG3", "CD28", "IL7R"))
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/t-reg_ftr_plot.png"), 
            plot = p)

    #Naive CD8 feature plot: Done
    p <- FeaturePlot(seurat.object, features = c("CD8A", "LEF1", "SELL", "CCR7", 
                                                "IL7R", "TCF7", "FOXP1", "BCL2",
                                                "CD28", "PTPN7"))
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/naive_cd8_ftr_plot.png"))

    #Memory CD8 T cell feature plot: Done
    p <- FeaturePlot(seurat.object, features = c("CD8A", "IL2RB", "KLRG1", "IL7R", #Effector memory marker 
                                                "TBX21", "PRDM1", "ID2", "STAT4", "NKG7", #Effector memory markers
                                                "SELL", "CCR7","EOMES", "TCF7", #Central memory markers
                                                "BCL6", "ID3", "STAT3")) #Central memory marker
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/memory_cd8_ftr_plot.png"), 
            plot = p)

    #MAIT cell feature plot: Done
    p <- FeaturePlot(seurat.object, features = c("CD8A", "IL18R1", "DPP4", "KLRB1", #Upregulated
                                                "NKG7", "ZBTB16", "SLC4A10", "ABCB1", "MAF", #Upregulated
                                                "CCR7", "TBX21", "SELL")) #Downregulated
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/MAIT_ftr_plot.png"), 
            plot = p)


    #gdT-cell feature plot
    p <- FeaturePlot(seurat.object, features = c("KLRD1", "TRDC", "ZAP70", "NKG7",
                                                "NCR1", "PRF1", "GZMB", "IL2RB", "CD247")) 
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/gdT_ftr_plot.png"), 
            plot = p)

    #Vd1 gdT-cell feature plot
    p <- FeaturePlot(seurat.object, features = c("CACHD1", "SELL", "IL7R", "CCR7", 
                                                "CX3CR1", "GZMA", "PRF1")) 
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/Vd1_gdT_ftr_plot.png"), 
            plot = p)

    #Vd2 gdT-cell feature plot
    p <- FeaturePlot(seurat.object, features = c("KLRD1", "GZMB", "PRF1", "IFNG", 
                                                "NKG7", "S100A4", "FCGR3A", "GNLY", 
                                                "CST7")) 
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/Vd2_gdT_ftr_plot.png"), 
            plot = p)

    #NK cell feature plot: Done
    p <- FeaturePlot(seurat.object, features = c("NCAM1", "ITGAM", "GZMB", "PRF1",
                                                "FCGR3A", "ITGB2", "KLRD1", "KLRC1",
                                                "KIT", "IL2RB", "SH2D1B"))
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/NK_ftr_plot.png"), 
            plot = p)

    #B cell feature plot: Done
    p <- FeaturePlot(seurat.object, features = c("CD19", "MS4A1", "CD79A", "PAX5",
                                                "IGHD", "IGHM", "TNFRSF13C", "BLK",
                                                "CD22"))
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/b_cell_ftr_plot.png"), 
            plot = p)
    
    # Precursor markers
    p <- FeaturePlot(seurat.object, features = c("CD34", "KIT", "CD38", "MKI67",
                                                 "CD33", "CD13"))
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/precursor_ftr_plot.png"), 
           plot = p)
}