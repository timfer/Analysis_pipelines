#Scope | 
# 1) Perform Seurat batch correction
# 2) Batch correction experiment-based
#Script-Type | Function
#Input | count matrices
#Output | seurat object
#Author | Joern Pezoldt
#Date | 30.05.2024

#####
#Parameter examples
#####
#species = "human"
#integrate_over = "CC"

#Threshing
#min_cells_percent = 0.005
#min_gene_number = 1000
#mito_cutoff = 15
#min_nUMI = 3000
##Seurat
#dims_use = 30
#resolution = 0.3

#common variables to regress are: "n_UMI", "percent.mt", "Phase"
#vars_to_regress <- c()

custom.palette <- c("1"="#104E8B", "2"="#FF3030", "3"="#228B22", "4"="#FF8C00",
                    "5"="#9400D3", "6"="#00FFFF", "7"="#FF1493", "8"="#7FFF00",
                    "9"="#FFD700", "10"="#FF00FF", "11"="#00FF7F", "12"="#FF0000",
                    "13"="#FFFF00", "14"="#D02090", "15"="#A020F0", "16"="#ADFF2F",
                    "17"="#483D8B", "18"="#00CED1", "19"="#00FF00", "20"="#FFA500",
                    "21"="#8B008B", "22"="#FF4500", "23"="#E066FF", "24"="#00BFFF",
                    "25"="#FF69B4", "26"="#CD0000", "27"="#54FF9F", "30"="#EE1289")

#Conditionals
# 1) If-clause for species for cell-cycle genes, mt, rpl, hsp, computation
seurat_integrate <- function(
  signature_cell_cycle_human_mouse, #signature cell cycle
  species,#species
  cond.2.rmv, #Conditions to remove
  matrix, #UMI matrix
  meta_data, #meta data
  #Variables
  integrate_over, # Batch over "CC" or "expID"
  #Seurat variables
  min_cells_percent,
  min_gene_number,
  mito_cutoff,
  dims_use,
  vars_to_regress, #common variables to regress are: "n_UMI", "percent.mt", "Phase"
  split.layers, #Variable for layer splitting and dowsntream integration
  integrate.layers,
  integration.type,
  #Path to save figure
  path) {
  #####
  #Load libraries
  #####
  require(Matrix)
  require(Seurat)
  #####
  #Setup species intel
  #####
  if(species == "human"){
    signature_cell_cycle <- subset(signature_cell_cycle_human_mouse,
                                   species == "human")
    mito_genes = "^MT-"
    ribo_genes = "^RPL"
    heatshock_genes = "^HSP"
  }
  if(species == "mouse"){
    signature_cell_cycle <- subset(signature_cell_cycle_human_mouse,
                                   species == "mouse")
    mito_genes = "^mt-"
    ribo_genes = "^Rpl"
    heatshock_genes = "^Hsp"
  }
  print("<<<Species intel set>>>")
  
  ######
  #Run Seurat
  ######
  matrix <- Matrix(matrix, sparse = TRUE)
  #Ensures that rownames are not considered duplicated
  rownames(matrix) <- make.unique(rownames(matrix))
  print("<<<Count matrix>>>")
  
  #Create Seurat objects
  ##Compute minimum number of cells expressing a gene
  min_cells <- ncol(matrix) * min_cells_percent
  seurat.object_v5 <- CreateSeuratObject(matrix, assay = "RNA", 
                                         meta.data = meta_data,
                                         min.cells = min_cells, 
                                         min.features = min_gene_number)
  #Remove unwanted conditions
  seurat.object_v5 <- subset(seurat.object_v5, 
                             cells = which(!(seurat.object_v5@meta.data$expID %in% cond.2.rmv)))
  print("Experiments kept: ")
  print(unique(seurat.object_v5@meta.data$expID))
  
  #Extract factor over which to integrate
  if(split.layers == "yes"){ 
    integrate_by <- seurat.object_v5[[integrate_over]]
    integrate_by_FCI <- rownames(integrate_by)
    integrate_by_type <- integrate_by[,1]
    names(integrate_by_type) <- integrate_by_FCI
    
    #Split data
    seurat.object_v5[["RNA"]] <- split(seurat.object_v5[["RNA"]], 
                                       f = integrate_by_type)
    print(paste("<<<Seurat layers split according to integration parameter ",
                integrate_over,">>>"))
  } else {
    print(paste("Count matrix not split according to ", integrate_over))
  }
  
  #Calculate QC genes
  seurat.object_v5[["percent.mt"]] <- PercentageFeatureSet(seurat.object_v5, 
                                                           pattern = mito_genes)
  seurat.object_v5[["percent.Rpl"]] <- PercentageFeatureSet(seurat.object_v5, 
                                                            pattern = ribo_genes)
  seurat.object_v5[["percent.Hsp"]] <- PercentageFeatureSet(seurat.object_v5, 
                                                            pattern = heatshock_genes)
  #Filter cells
  print(paste("All :",nrow(seurat.object_v5@meta.data)))
  seurat.object_v5 <- subset(seurat.object_v5, subset = nFeature_RNA  > min_gene_number)
  print(paste("After Feature Thresh :",nrow(seurat.object_v5@meta.data)))
  seurat.object_v5 <- subset(seurat.object_v5, subset = percent.mt < mito_cutoff)
  print(paste("After Mito Thresh :",nrow(seurat.object_v5@meta.data)))
  print("<<<Seurat object threshed for cutoffs>>>")
  
  seurat.object_v5 <- NormalizeData(seurat.object_v5)
  #DefaultLayer(seurat.object_v5[["RNA"]]) <- 'data'
  
  #Find variable genes, scale and PCA
  seurat.object_v5 <- FindVariableFeatures(seurat.object_v5, verbose = FALSE)
  print("<<<Variable features found>>>")
  
  print("<<<Scaling data>>>")
  seurat.object_v5 <- ScaleData(seurat.object_v5,
                                assay = "RNA",
                                model.use = "linear",
                                features = rownames(seurat.object_v5@assays$RNA),
                                vars.to.regress = vars_to_regress,
                                verbose = TRUE)
  print("<<<Data scaled>>>")
  
  seurat.object_v5 <- RunPCA(seurat.object_v5, npcs = dims_use, verbose = TRUE)
  integration.name <- "pca"
  
  # Ensure the output directory exists
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  
  #Generate elbow plot for dimension determination
  p <- ElbowPlot(seurat.object_v5, ndims = dims_use) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_line(color = "grey90")
    )
  print(p)
  ggsave(filename = paste0(path, "/elbow_plot.png"), plot = p)
  print("<<<FVB, Scale, Regression, PCA done>>>")
  
  if (!dir.exists(paste0(path, "/PCA_hmap"))) {
    dir.create(paste0(path, "/PCA_hmap"), recursive = TRUE)
  }
  
  # Loop over blocks of 10 PCs
  for (i in 1:14) {
    dims_to_plot <- (i * 10):((i + 1) * 10)
    fname <- sprintf("dim_%d-%d.png", dims_to_plot[1], dims_to_plot[length(dims_to_plot)])
    fpath <- file.path(paste0(path, "/PCA_hmap/", fname))
    
    # Open a PNG device (8x12 inches at 300 dpi)
    png(
      filename = fpath,
      width    = 8,
      height   = 12,
      units    = "in",
      res      = 300
    )
    
    # This draws the heatmap into the file
    DimHeatmap(
      object   = seurat.object_v5,
      dims     = dims_to_plot,
      cells    = 1000,
      balanced = TRUE
    )
    
    # Close the device, finalizing the PNG
    dev.off()
  }
  
  # Score and plot JackStraw
  print("<<<Scoring and plotting JackStraw>>>")
  # max_dims <- min(dims_use, length(seurat.object_v5@reductions$pca))

  seurat.object_v5 <- JackStraw(seurat.object_v5, reduction = "pca", dims=dims_use,
                               num.replicate = 100, verbose = TRUE)
  seurat.object_v5 <- ScoreJackStraw(seurat.object_v5, reduction = "pca",
                                     dims = 1:150)
  # Generate and save JackStrawPlot
  p <- JackStrawPlot(seurat.object_v5, dims = 1:150) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_line(color = "grey90")
    )
  print(p)
  ggsave(filename = paste0(path, "/jackstraw_plot.png"), plot = p)
  
  p <- JackStrawPlot(seurat.object_v5, dims = 100:150) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_line(color = "grey90")
    )
  print(p)
  ggsave(filename = paste0(path, "/jackstraw_plot_last_dims.png"), plot = p)
  
  #Compute Integrations
  # Ensure k.weight is smaller than the minimum number of cells
  k.weight <- floor(min(100, min_cells))
  cat("k.weight set to ", k.weight, "\n")
  
  if(integrate.layers == "yes"){
    if(integration.type == "cca"){
      integration.name <- "integrated.cca"
      seurat.object_v5 <- IntegrateLayers(
        object = seurat.object_v5, method = CCAIntegration,
        orig.reduction = "pca", new.reduction = integration.name,
        k.weight = k.weight,
        verbose = FALSE)
      print(paste("<<<Multiple layers integrated according to CCAIntegration>>>"))
    }else if(integration.type == "rpca"){
      integration.name <- "integrated.rpca"
      seurat.object_v5 <- IntegrateLayers(
        object = seurat.object_v5, method = RPCAIntegration,
        orig.reduction = "pca", new.reduction = integration.name,
        k.weight = k.weight,
        verbose = FALSE)
      print(paste("<<<Multiple layers integrated according to RPCAIntegration>>>"))
    }else if(integration.type == "harmony"){
      integration.name <- "integrated.harmony"
      # 1. Run Harmony on the PCA embeddings
      seurat.object_v5 <- RunHarmony(
        object       = seurat.object_v5,
        group.by.vars = integrate_over,       # or your batch variable
        reduction     = "pca",
        assay         = DefaultAssay(seurat.object_v5),
        project.dim   = FALSE,
        reduction.save = integration.name
      )
      print(paste("<<<Multiple layers integrated according to Harmony>>>"))
    }
  }
  
  #Find neighbors
  seurat.object_v5 <- FindNeighbors(seurat.object_v5, reduction = integration.name, 
                                    dims = 1:dims_use, verbose = TRUE)
  print("<<<Neighbors found>>>")
  return(list(seurat.object_v5, dims_use, integration.name))
}

find_cluster_resolution <- function(seurat.object_v5, 
                                    dims_use, resolutions, min.pct, logfc.thresh,
                                    integration.name, algorithm, path){
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  
  print(paste("Algorithm used: ", algorithm))
  if(algorithm == "Leiden"){
    algo = 4
  } else if(algorithm == "Louvain") {
    algo = 1
  }

  for(res in resolutions){
    reduction_name <- paste0(integration.name, "umap.res", res)
    # Run clustering with different resolutions
    seurat.object_v5 <- FindClusters(seurat.object_v5, resolution = res, algorithm = algo,
                                     cluster.name = "cca_clusters", verbose = TRUE)
    #Run UMAP with current resolution
    seurat.object_v5 <- RunUMAP(seurat.object_v5, reduction = integration.name,
                                dims = 1:dims_use, reduction.name = reduction_name, 
                                verbose = TRUE)
    p <- DimPlot(seurat.object_v5, reduction = reduction_name, label = TRUE, pt.size = 0.3,
                 cols = custom.palette) +
      labs(title = paste("UMAP clustering at resolution", res)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    print(p)
    ggsave(filename = paste0(path, "/", gsub("\\.", "_", reduction_name), ".png"), 
           plot = p, height = 8, width = 10)
    
    seurat.object_v5 <- JoinLayers(seurat.object_v5)
    markers.all <- FindAllMarkers(seurat.object_v5, only.pos = FALSE,
                                  min.pct = min.pct, logfc.threshold = logfc.threshold)

    markers.all <- markers.all %>% arrange(cluster, desc(avg_log2FC))
    head(markers.all)

    #Select top 250 DEGs
    markers.all %>% group_by(cluster) %>%
      top_n(n = 7, wt = avg_log2FC) -> top7

    #Plot heatmap
    ordered_palette <- custom.palette[order(as.integer(names(custom.palette)))]
    p <- DoHeatmap(seurat.object_v5, features = top7$gene,
                   size = 3, angle = 90, group.colors = ordered_palette) +
      scale_fill_gradient2( low = rev(c('#0000CD','#000080', '#00008B')),
                            mid = "white", high = rev(c('#FFD700','#EEC900','#CDAD00')),
                            midpoint = 0, guide = "colourbar", aesthetics = "fill") +
      NoLegend()
    print(p)
    ggsave(filename = paste0(path, "/", gsub("\\.", "_", paste0(integration.name, "DEG_res", res)), ".png"),
           plot = p, width = 10, height = 18)
  }
}

cluster_UMAP_seurat <- function(seurat.object, resolution, dims, 
                                integration.name, reduction.name, algorithm, path){
  print(paste("Algorithm used: ", algorithm))
  if(algorithm == "Leiden"){
    algo = 4
  } else if(algorithm == "Louvain") {
    algo = 1
  }
  #Determine UMAP clusters
  seurat.object <- FindClusters(seurat.object, resolution = resolution, algorithm = algo,
                                cluster.name = paste0(integration.name, "_clusters"), verbose = TRUE)
  #Run UMAP with current resolution
  seurat.object <- RunUMAP(seurat.object, reduction = integration.name,
                           dims = 1:dims, reduction.name = reduction.name, 
                           verbose = TRUE, return.model = TRUE)
  print("<<<UMAP coordinates determined>>>")
  
  return(seurat.object)
}
