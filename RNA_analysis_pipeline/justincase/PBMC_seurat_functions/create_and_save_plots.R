# Set custom color palette here

generate_QC_plots <- function(seurat.object, res, reduction, target_expID, path) {
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
  ggsave(filename = paste0(path, "/nftr_RNA_vln_expID.pdf"), plot = p,
         width = 8, height = 8, device = 'pdf', units = "in", useDingbats = FALSE)
  
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
  ggsave(filename = paste0(path, "/ncount_RNA_vln_expID.pdf"), plot = p,
         width = 8, height = 8, device = 'pdf', units = "in", useDingbats = FALSE)
  
  p <- VlnPlot(seurat.object, features = c("percent.mt"),
           group.by = "expID", pt.size = 0) + 
    labs(title = NULL) +
    xlab("Experimental ID (batch)") + ylab("Mitochondrial gene expression (per cent (%))") + 
    plot_annotation(title = "Fraction of mitochondrial genes expressed per experimental batch") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    NoLegend() +
    geom_boxplot(width = 0.1, fill = NA, color = "black", position = position_dodge(width = 0.9))
  print(p)
  ggsave(filename = paste0(path, "/pct_mt_vln_expID.pdf"), plot = p,
         width = 8, height = 8, device = 'pdf', units = "in", useDingbats = FALSE)
  
  p <- VlnPlot(seurat.object, features = c("nFeature_RNA"),
               group.by = "seurat_clusters", pt.size = 0) + 
    labs(title = NULL) +
    xlab(paste0("UMAP cluster (resolution ", res,")")) + ylab("nFeature_RNA (counts)") + 
    plot_annotation(title = "Feature RNA counts per integrated low-dimensional cluster") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    NoLegend() +
    geom_boxplot(width = 0.1, fill = NA, color = "black", position = position_dodge(width = 0.9))
  print(p)
  ggsave(filename = paste0(path, "/nftr_RNA_vln_cluster.pdf"), plot = p,
         width = 8, height = 8, device = 'pdf', units = "in", useDingbats = FALSE)
  
  p <- VlnPlot(seurat.object, features = c("nCount_RNA"),
               group.by = "seurat_clusters", pt.size = 0) + 
    labs(title = NULL) +
    xlab(paste0("UMAP cluster (resolution ", res,")")) + ylab("nCount_RNA (counts)") + 
    plot_annotation(title = "RNA counts per integrated low-dimensional cluster") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    NoLegend() +
    geom_boxplot(width = 0.1, fill = NA, color = "black", position = position_dodge(width = 0.9))
  print(p)
  ggsave(filename = paste0(path, "/ncount_RNA_vln_cluster.pdf"), plot = p,
         width = 8, height = 8, device = 'pdf', units = "in", useDingbats = FALSE)
  
  p <- VlnPlot(seurat.object, features = c("percent.mt"),
               group.by = "seurat_clusters", pt.size = 0) + 
    labs(title = NULL) +
    xlab(paste0("UMAP cluster (resolution ", res,")")) + ylab("Mitochondrial gene expression (per cent (%))") + 
    plot_annotation(title = "Feature RNA counts per integrated low-dimensional cluster") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    NoLegend() +
    geom_boxplot(width = 0.1, fill = NA, color = "black", position = position_dodge(width = 0.9))
  print(p)
  ggsave(filename = paste0(path, "/pct_mt_vln_cluster.pdf"), plot = p,
         width = 8, height = 8, device = 'pdf', units = "in", useDingbats = FALSE)
  
  p1 <- FeatureScatter(seurat.object, feature1 = "n_mapped_reads", feature2 = "n_UMI",
                       group.by = "expID") +
    scale_x_continuous(labels = function(x) x / 1e4) + xlab(expression("n_mapped_reads (×10"^4*")")) +
    scale_y_continuous(labels = function(x) x / 1e2) + ylab(expression("n_UMI (×10"^2*")"))
  p2 <- FeatureScatter(seurat.object, feature1 = "n_mapped_reads", feature2 = "nFeature_RNA",
                       group.by = "expID") +
    scale_x_continuous(labels = function(x) x / 1e4) + xlab(expression("n_mapped_reads (×10"^4*")"))
  p <- p1 + p2 + plot_layout(guides = 'collect')
  print(p)
  ggsave(filename = paste0(path, "/reads_ftr_scatter.pdf"), plot = p,
         width = 8, height = 8, device = 'pdf', units = "in", useDingbats = FALSE)
  
  p1 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA",
                      feature2 = "percent.mt", group.by = "expID")
  p2 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA",
                      feature2 = "nFeature_RNA", group.by = "expID")
  p <- p1 + p2 + plot_layout(guides = 'collect')
  print(p)
  ggsave(filename = paste0(path, "/var_ftr_scatter.pdf"), plot = p,
         width = 8, height = 8, device = 'pdf', units = "in", useDingbats = FALSE)
  
  #Plot regression variables on feature scatter plot
  ftrs = c("nCount_RNA","nFeature_RNA", "percent.mt", "percent.Hsp")
  p <- FeaturePlot(seurat.object, features = ftrs, reduction = reduction, 
                  pt.size = 0.3)
  print(p)
  ggsave(filename = paste0(path, "/ftr_plot_regression_vars.pdf"), plot = p,
         width = 16, height = 16, device = 'pdf', units = "in", useDingbats = FALSE)
  
  # #Plot a single target experiment vs the rest
  # if(any(target_expID != "")){
  #   plots_list <- list()
  #   for(expID_i in target_expID){
  #     target_subset <- subset(seurat.object, expID == expID_i)
  #     p_target <- DimPlot(target_subset, reduction = reduction,
  #                         label = FALSE, pt.size = 0.3) + #, cols = cluster_colors) +
  #       NoLegend() + ggtitle(paste("Experiment: ", expID_i))
  #     plots_list[[length(plots_list)+1]] <- p_target
  #     }
  #   others_subset <- subset(seurat.object, cells = which(!(seurat.object@meta.data$expID %in% target_expID)))
  #   p_others <- DimPlot(others_subset, reduction = reduction,
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
  # ggsave(filename = paste0(path, "/UMAP_target_exp_v_rest.pdf"), plot = p,
  #          width = width, height = 8, limitsize = FALSE)
  # } else {
  #   message("target_expID is empty, skipping the plot.")
  # }
}

generate_bar_plots <- function(seurat.object, path) {
  # Prepare data by grouping based on Seurat cluster, 'orig.ident' (experimental ID), and 'CC'
  cluster_counts <- seurat.object@meta.data %>%
    group_by(seurat_clusters, orig.ident, CC) %>%
    tally() %>% ungroup()
  
  # Absolute Bar Plot
  p <- ggplot(cluster_counts, aes(x = seurat_clusters, y = n, fill = CC)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    # scale_fill_manual(values = custom.palette) +
    facet_wrap(~ orig.ident) +
    labs(title = "Absolute Number of Cells by Cluster", 
         y = "Cell Count (absolute)", 
         x = "Seurat Cluster") # +
    # theme_minimal() +
    # theme(panel.background = element_rect(fill = "white", color = NA))
  
  print(p)
  ggsave(filename = paste0(path, "/absolute_distribution_bar_plot.pdf"), plot = p,
         width = 12, height = 8, device = 'pdf', units = "in", useDingbats = FALSE)
  
  cluster_counts <- cluster_counts %>%
    group_by(orig.ident, CC) %>%
    mutate(relative_count = n / sum(n)) %>%
    ungroup()
  
  # Relative Bar Plot
  p <- ggplot(cluster_counts, aes(x = seurat_clusters, y = relative_count, fill = CC)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    # scale_fill_manual(values = custom.palette) +
    facet_wrap(~ orig.ident) +
    labs(title = "Relative Proportion of Cells by Cluster", 
         y = "Relative Cell Count (%)", 
         x = "Seurat Cluster") # +
    # theme_minimal() +
    # theme(panel.background = element_rect(fill = "white", color = NA))
  
  print(p)
  ggsave(filename = paste0(path, "/relative_distribution_bar_plot.pdf"), plot = p,
         width = 12, height = 8, device = 'pdf', units = "in", useDingbats = FALSE)
}

plot_specific_ftr <- function(seurat.object, path){
    ftr_plot_path = path
    if (!dir.exists(ftr_plot_path)) {
            dir.create(ftr_plot_path, recursive = TRUE)
    }

    #Monocyte feature plot: Done
    p <- FeaturePlot(seurat.object, features = c("HLA.DRA", "CD14", "CSF3R", "CCR2", #Classical monocytes
                                                "LYZ", "S100A8", "S100A9", "ITGAM",#Classical monocytes
                                                "FCGR3A", "CX3CR1", "CD86", "ITGB2", "ITGAL"),
                     ncol = 4) #Non-classical monocytes
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/monocyte_ftr_plot.pdf"), 
            plot = p, width = 16, height = 20, 
           device = 'pdf', units = "in", useDingbats = FALSE)

    #Naive CD4 T-cell feature plot: Done
    p <- FeaturePlot(seurat.object, features = c("CD4", "LEF1", "TCF7", "CCR7", 
                                                 "SELL", "LTB", "IL7R", "CD28", 
                                                 "FOXP1", "CD3D", "CD69", "BCL2"),
                     ncol = 4) 
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/naive_cd4_ftr_plot.pdf"), 
            plot = p, width = 16, height = 12, 
           device = 'pdf', units = "in", useDingbats = FALSE)

    #Central memory CD4+ T cells
    p <- FeaturePlot(seurat.object, features = c("CD4", "IL2RA", "CD44", "SELL"),
                     ncol = 2)
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/central_mem_cd4_ftr_plot.pdf"), 
            plot = p, width = 8, height = 8, 
           device = 'pdf', units = "in", useDingbats = FALSE)

    #Th2 effector memory CD4 feature plot: Done
    p <- FeaturePlot(seurat.object, features = c("GATA3", "STAT6", "STAT5A", 
                                                 "CCR4", "RELA", "CD40LG",
                                                 "ICOS", "TNFRSF4"),
                     ncol = 3)
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/memory_cd4_Th2_ftr_plot.pdf"), 
           plot = p, width = 12, height = 12, 
           device = 'pdf', units = "in", useDingbats = FALSE)

    #T-reg CD4 feature plot: Done
    p <- FeaturePlot(seurat.object, features = c("IL2RA", "CTLA4", 
                                                "IKZF2", "TGFB1", "ITGB1",
                                                "LAG3", "CD28", "IL7R"),
                     ncol = 3)
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/t-reg_ftr_plot.pdf"), 
            plot = p, width = 12, height = 12, device = 'pdf', units = "in", useDingbats = FALSE)

    #Naive CD8 feature plot: Done
    p <- FeaturePlot(seurat.object, features = c("CD8A", "LEF1", "SELL", 
                                                 "CCR7", "IL7R", "TCF7",
                                                 "FOXP1", "BCL2", "CD28",
                                                 "PTPN7"),
                     ncol = 3)
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/naive_cd8_ftr_plot.pdf"),
           width = 12, height = 16, device = 'pdf', units = "in", useDingbats = FALSE)

    #Memory CD8 T cell feature plot: Done
    p <- FeaturePlot(seurat.object, features = c("CD8A", "IL2RB", "KLRG1", "IL7R", #Effector memory marker 
                                                "TBX21", "PRDM1", "ID2", "STAT4", 
                                                "NKG7", #Effector memory markers
                                                "SELL", "CCR7","EOMES", 
                                                "TCF7", #Central memory markers
                                                "BCL6", "ID3", "STAT3"),#Central memory marker
                     ncol = 4)
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/memory_cd8_ftr_plot.pdf"), 
            plot = p, width = 16, height = 16, device = 'pdf', units = "in", useDingbats = FALSE)

    #MAIT cell feature plot: Done
    p <- FeaturePlot(seurat.object, features = c("CD8A", "IL18R1", "DPP4", "KLRB1", #Upregulated
                                                "NKG7", "ZBTB16", "SLC4A10", "ABCB1", "MAF", #Upregulated
                                                "CCR7", "TBX21", "SELL"), #Downregulated
                     ncol = 4) 
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/MAIT_ftr_plot.pdf"), 
            plot = p, width = 16, height = 12, device = 'pdf', units = "in", useDingbats = FALSE)


    #gdT-cell feature plot
    p <- FeaturePlot(seurat.object, features = c("KLRD1", "TRDC", "ZAP70", 
                                                 "NKG7", "NCR1", "PRF1", 
                                                 "GZMB", "IL2RB", "CD247"),
                     ncol = 3) 
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/gdT_ftr_plot.pdf"), 
            plot = p, width = 12, height = 12, device = 'pdf', units = "in", useDingbats = FALSE)

    #Vd1 gdT-cell feature plot
    p <- FeaturePlot(seurat.object, features = c("CACHD1", "SELL", "IL7R", 
                                                 "CCR7", "CX3CR1", "GZMA", 
                                                 "PRF1"), ncol = 3) 
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/Vd1_gdT_ftr_plot.pdf"), 
            plot = p, width = 12, height = 12, device = 'pdf', units = "in", useDingbats = FALSE)

    #Vd2 gdT-cell feature plot
    p <- FeaturePlot(seurat.object, features = c("KLRD1", "GZMB", "PRF1", 
                                                 "IFNG", "NKG7", "S100A4", 
                                                 "FCGR3A", "GNLY", "CST7"),
                     ncol = 3) 
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/Vd2_gdT_ftr_plot.pdf"), 
            plot = p, width = 12, height = 12, device = 'pdf', units = "in", useDingbats = FALSE)

    #NK cell feature plot: Done
    p <- FeaturePlot(seurat.object, features = c("NCAM1", "ITGAM", "GZMB", "PRF1",
                                                "FCGR3A", "ITGB2", "KLRD1", "KLRC1",
                                                "KIT", "IL2RB", "SH2D1B"), ncol = 4)
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/NK_ftr_plot.pdf"), 
            plot = p, width = 16, height = 12, device = 'pdf', units = "in", useDingbats = FALSE)

    #B cell feature plot: Done
    p <- FeaturePlot(seurat.object, features = c("CD19", "MS4A1", "CD79A", 
                                                 "PAX5", "IGHD", "IGHM",
                                                 "TNFRSF13C", "BLK", "CD22"), ncol = 3)
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/b_cell_ftr_plot.pdf"), 
            plot = p, width = 12, height = 12, device = 'pdf', units = "in", useDingbats = FALSE)
    
    # Precursor markers
    p <- FeaturePlot(seurat.object, features = c("CD34", "KIT", "CD38", 
                                                 "MKI67", "CD33", "ANPEP", 
                                                 "MPO", "BCL2"), ncol = 3)
    print(p)
    ggsave(filename = paste0(ftr_plot_path, "/leukemia_ftr_plot.pdf"), 
           plot = p, width = 12, height = 12, device = 'pdf', units = "in", useDingbats = FALSE)
}