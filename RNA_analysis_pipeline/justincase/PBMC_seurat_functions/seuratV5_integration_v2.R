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
  reduction,
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
    print(paste("<<<Seurat layers split according to integration parameter",integrate_over,">>>"))
  } else {
    print(paste("Count matrix not split according to",integrate_over))
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
  
  # Ensure the output directory exists
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  
  # Generate elbow plot for dimension determination
  p <- ElbowPlot(seurat.object_v5, ndims = dims_use) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_line(color = "grey90")
    )
  print(p)
  ggsave(filename = paste0(path, "/elbow_plot.pdf"), plot = p,
         device = 'pdf', units = "in", useDingbats = FALSE)
  print("<<<FVB, Scale, Regression, PCA done>>>")

  if (!dir.exists(paste0(path, "/PCA_hmap"))) {
    dir.create(paste0(path, "/PCA_hmap"), recursive = TRUE)
  }
  
  # Loop over blocks of 10 PCs
  block_size <- 10
  start_positions <- seq(1, dims_use, by = block_size)
  for (i in start_positions) {
    # Define the range for this iteration
    end_pos = min(i + block_size - 1, dims_use)
    dims_to_plot = i:end_pos
    
    fname <- sprintf("dim_%d-%d.pdf", dims_to_plot[1], dims_to_plot[length(dims_to_plot)])
    hmappath <- file.path(paste0(path, "/PCA_hmap/", fname))
    
    # Open a PNG device (8x12 inches at 300 dpi)
    png(
      filename = hmappath,
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
  print("<<<Scoring and plotting JackStraw>>>")
  seurat.object_v5 <- JackStraw(seurat.object_v5, reduction = "pca", dims=dims_use,
                                num.replicate = 100, verbose = TRUE)
  seurat.object_v5 <- ScoreJackStraw(seurat.object_v5, reduction = "pca",
                                     dims = 1:dims_use)
  
  if (!dir.exists(paste0(path, "/jackstraw_plots"))) {
    dir.create(paste0(path, "/jackstraw_plots"), recursive = TRUE)
  }
  for (i in start_positions) {
    # Define the range for this iteration
    end_pos = min(i + block_size - 1, dims_use)
    dims_to_plot = i:end_pos
    
    fname <- sprintf("dim_%d-%d.pdf", dims_to_plot[1], dims_to_plot[length(dims_to_plot)])
    jspath <- file.path(paste0(path, "/jackstraw_plots/", fname))
    p <- JackStrawPlot(seurat.object_v5, dims = dims_to_plot) +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey90")
      )
    print(p)
    ggsave(filename = jspath, plot = p, width = 8, height = 8, 
           device = 'pdf', units = "in", useDingbats = FALSE)
    
  }

  #Compute Integrations
  # Ensure k.weight is smaller than the minimum number of cells
  k.weight <- floor(min(100, min_cells))
  cat("k.weight set to ", k.weight, "\n")
  
  if(integrate.layers == "yes"){
    if(reduction == "cca"){
      seurat.object_v5 <- IntegrateLayers(
        object = seurat.object_v5, method = CCAIntegration,
        orig.reduction = "pca", new.reduction = reduction,
        k.weight = k.weight,
        verbose = FALSE)
      print(paste("<<<Multiple layers integrated according to CCAIntegration>>>"))
    }else if(reduction == "rpca"){
      seurat.object_v5 <- IntegrateLayers(
        object = seurat.object_v5, method = RPCAIntegration,
        orig.reduction = "pca", new.reduction = reduction,
        k.weight = k.weight,
        verbose = FALSE)
      print(paste("<<<Multiple layers integrated according to RPCAIntegration>>>"))
    }else if(reduction == "harmony"){
      # 1. Run Harmony on the PCA embeddings
      seurat.object_v5 <- RunHarmony(
        object       = seurat.object_v5,
        group.by.vars = integrate_over,       # or your batch variable
        reduction     = "pca",
        assay         = DefaultAssay(seurat.object_v5),
        project.dim   = FALSE,
        reduction.save = reduction
      )
      print(paste("<<<Multiple layers integrated according to Harmony>>>"))
    }
  }
  return(list(seurat.object_v5, dims_use))
}

cluster_UMAP_seurat <- function(seurat.object, resolution, dims, 
                                reduction, algorithm, path){
  print(paste("Algorithm used: ", algorithm))
  if(algorithm == "Leiden"){
    algo = 4
  } else if(algorithm == "Louvain") {
    algo = 1
  }
  #Determine UMAP clusters
  seurat.object <- FindClusters(seurat.object, resolution = resolution, algorithm = algo,
                                cluster.name = paste0(reduction, "_clusters"), verbose = TRUE)
  #Run UMAP with current resolution
  seurat.object <- RunUMAP(seurat.object, reduction = reduction,
                           dims = 1:dims, verbose = TRUE, return.model = TRUE)
  print(paste0("<<<UMAP coordinates determined at resolution:", resolution , ">>>"))
  
  if (!dir.exists(paste0(path, "/split_umaps"))) {
    dir.create(paste0(path, "/split_umaps"), recursive = TRUE)
  }
  ######
  # Plot the experimental overlap/splits with the final cluster resolutions
  ######
  #Plot UMAP with overlayed experiments to check for batch effects
  # colors <- brewer.pal(length(unique(seurat.object$expID)), "Set2") #Set colors
  p <- DimPlot(seurat.object, reduction = "umap", pt.size = 0.5,
               group.by = "expID")
    print(p)
  ggsave(filename = paste0(path, "/split_umaps", "/UMAP_exp_overlap.pdf"), plot = p,
         width = 16, height = 8, device = 'pdf', units = "in", useDingbats = FALSE)
  
  # Replot split by experiment now that dimensions and resolutions have been picked
  p <- DimPlot(seurat.object, reduction = "umap", label = FALSE,
               pt.size = 0.5, split.by = "expID", ncol = 4) + NoLegend()
  print(p)
  ggsave(
    filename = paste0(path, "/split_umaps", "/UMAP_exp_split.pdf"), 
    plot = p,
    width = 16, 
    height = ceiling(length(unique(seurat.object_qc@meta.data$expID)) / 4) * 4, 
    limitsize = FALSE, device = 'pdf', units = "in", useDingbats = FALSE
  )
  
  # Plot the UMAP
  
  p <- DimPlot(seurat.object, label = TRUE) #, cols = umap_colors)
  print(p)
  ggsave(filename = paste0(path, "/UMAP_unnamed.pdf"), plot = p,
         height = 8, width = 10)
  
  #####
  #Feature plots & UMAP
  #####
  #Plot cell type discriminating genes on feature plot
  ftrs <- c("CD3E", "CD8A", "CD4", "CD19", 
           "CD14", "NCAM1", "NKG7", "TNF", 
           "IL6", "IL1B", "GZMB", "IL2RA")
  p <- FeaturePlot(seurat.object, features = ftrs, reduction = "umap",
                   pt.size = 0.3, ncol = 4)
  print(p)
  ggsave(filename = paste0(path, "/std_ftr_genes.pdf"), plot = p, width = 32, 
         height = 24, device = 'pdf', units = "in", useDingbats = FALSE)
  
   
  aml_ftrs <- c("NPM1", "KMT2A", "FLT3", "DNMT3A",
                "IDH1", "IDH2", "CEBPA", "SRSF2", 
                "TET2", "NRAS", "KRAS", "CD34", 
                "CD38", "KIT", "MKI67", "BCL2", 
                "MPO"
                )
  p <- FeaturePlot(seurat.object, features = aml_ftrs, reduction = "umap",
                   pt.size = 0.3, ncol = 4)
  print(p)
  ggsave(filename = paste0(path, "/aml_ftr_genes.pdf"), plot = p, width = 32, height = 40)
  
  mt_ff_ftrs <- c("MFN1", "MFN2", "OPA1",
                  "DNM1L", "MFF", "MIEF2", "MIEF1")
  p <- FeaturePlot(seurat.object, features = aml_ftrs, reduction = "umap",
                   pt.size = 0.3, ncol = 4)
  print(p)
  ggsave(filename = paste0(path, "/mt_ff_ftr_genes.pdf"), plot = p, width = 24, 
         height = 24, device = 'pdf', units = "in", useDingbats = FALSE)
  
  return(seurat.object)
}
