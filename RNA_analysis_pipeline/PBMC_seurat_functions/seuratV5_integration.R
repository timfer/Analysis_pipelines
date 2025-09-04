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
  #signature cell cycle
  signature_cell_cycle_human_mouse,
  #species
  species,
  cond.2.rmv, #Conditions to remove
  
  #UMI matrix
  matrix,
  #meta data
  meta_data,
  #Variables
  integrate_over, # Batch over "CC" or "expID"
  
  #Seurat
  min_cells_percent,
  min_gene_number,
  mito_cutoff,
  dims_use,
  resolution,
  #common variables to regress are: "n_UMI", "percent.mt", "Phase"
  vars_to_regress,
  join_layers,
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
    signature_cell_cycle <- subset(signature_cell_cycle_human_mouse, species == "human")
    mito_genes = "^MT-"
    ribo_genes = "^RPL"
    heatshock_genes = "^HSP"
  }
  if(species == "mouse"){
    signature_cell_cycle <- subset(signature_cell_cycle_human_mouse, species == "mouse")
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
  if(!is.null(cond.2.rmv)){
    seurat.object_v5 <- subset(seurat.object_v5, 
                               cells = which(!(seurat.object_v5@meta.data$expID %in% cond.2.rmv)))
    print("Experiments kept: ")
    print(unique(seurat.object_v5@meta.data$expID))
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
  
  #Cell cycle scoring
  s.genes <- subset(signature_cell_cycle, cycle_state == "S")$GeneSymbol
  g2m.genes <- subset(signature_cell_cycle, cycle_state == "G2/M")$GeneSymbol
  #s.genes <- toupper(s.genes)
  #g2m.genes <- toupper(g2m.genes)
  
  seurat.object_v5 <- NormalizeData(seurat.object_v5)
  DefaultLayer(seurat.object_v5[["RNA"]]) <- 'data'
  
  seurat.object_v5 <- CellCycleScoring(seurat.object_v5,
                                       s.features = s.genes,
                                       g2m.features = g2m.genes)
  print("<<<Cell cycle score computed>>>")
  
  #Extract factor over which to integrate
  if(!is.null(integrate_over)){
    integrate_by <- seurat.object_v5[[integrate_over]]
    integrate_by_FCI <- rownames(integrate_by)
    integrate_by_type <- integrate_by[,1]
    names(integrate_by_type) <- integrate_by_FCI
    
    #Split data
    seurat.object_v5[["RNA"]] <- split(seurat.object_v5[["RNA"]], 
                                       f = integrate_by_type)
    print("<<<Seurat layers according to integration parameter set>>>")
  }

  #Find variable genes, scale and PCA
  seurat.object_v5 <- FindVariableFeatures(seurat.object_v5, verbose = FALSE)
  
  seurat.object_v5 <- ScaleData(seurat.object_v5,
                                assay = "RNA",
                                model.use = "linear",
                                features = rownames(seurat.object_v5@assays$RNA),
                                vars.to.regress = vars_to_regress,
                                verbose = TRUE)

  seurat.object_v5 <- RunPCA(seurat.object_v5, npcs = dims_use, verbose = FALSE)
  
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
  
  stdevs <- seurat.object_v5@reductions$pca@stdev
  min_stdev <- min(stdevs)
  thresh <- min_stdev * 1.05 # Make a st_dev threshold @ 5% deviation
  #Return the dimension which fulfills the above threshold criterium, assign value to parameters
  dims_use <- which(stdevs <= thresh)[1]
  cat("PCA dimensions set: ", dims_use, "\n")
  print("<<<FVB, Scale, Regression, PCA done>>>")
  
  #Compute Integrations
  # Ensure k.weight is smaller than the minimum number of cells
  k.weight <- floor(min(100, min_cells))
  cat("k.weight set to ", k.weight, "\n")
  
  #Check number of count matrices in Layers
  layer_names <- names(seurat.object_v5@assays$RNA@layers)
  count_layers <- grep("count", layer_names, value = TRUE)
  
  if(length(count_layers) > 1){
    seurat.object_v5 <- IntegrateLayers(
      object = seurat.object_v5, method = CCAIntegration,
      orig.reduction = "pca", new.reduction = "integrated.cca",
      k.weight = k.weight,
      verbose = FALSE)
    print("<<<Multiple layers integrated>>>")
  }
  
  #Rename the "pca" reduction as "integrated.cca" if there is only 1 layer
  if(length(seurat.object_v5@reductions) == 1 && "pca" %in% names(seurat.object_v5@reductions)){
    # Rename the reduction
    seurat.object_v5@reductions[["integrated.cca"]] <- seurat.object_v5@reductions[["pca"]]
    seurat.object_v5@reductions[["pca"]] <- NULL
    
    print("Reduction renamed from 'pca' to 'integrated.cca'.")
  } else {
    print("Renaming not performed. Either there are multiple reductions or 'pca' does not exist.")
  }
  
  if(join_layers == "yes"){
    seurat.object_v5 <- JoinLayers(seurat.object_v5)
    print("<<<Layers joined>>>")
  }
  
  #Find neighbors
  seurat.object_v5 <- FindNeighbors(seurat.object_v5, reduction = "integrated.cca", 
                                    dims = 1:dims_use, verbose = TRUE)
  print("<<<Neighbors found>>>")
  return(list(seurat.object_v5, dims_use))
}

find_cluster_resolution <- function(seurat.object_v5, dims_use, 
                                    resolutions, path){
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  
  for(res in resolutions){
    reduction_name <- paste0("umap.cca.res", res)
    # Run clustering with different resolutions
    seurat.object_v5 <- FindClusters(seurat.object_v5, resolution = res,
                                     cluster.name = "cca_clusters", verbose = TRUE)
    #Run UMAP with current resolution
    seurat.object_v5 <- RunUMAP(seurat.object_v5, reduction = "integrated.cca",
                                dims = 1:dims_use, reduction.name = reduction_name, 
                                verbose = TRUE)
    p <- DimPlot(seurat.object_v5, reduction = reduction_name, label = TRUE, pt.size = 0.3)
    print(p)
    ggsave(filename = paste0(path, "/", gsub("\\.", "_", reduction_name), ".png"), plot = p)
  }
}

cluster_UMAP_seurat <- function(seurat.object, resolution, dims, reduction.name){
  #Determine UMAP clusters
  seurat.object <- FindClusters(seurat.object, resolution = resolution,
                                cluster.name = "cca_clusters", verbose = TRUE)
  #Run UMAP with current resolution
  seurat.object <- RunUMAP(seurat.object, reduction = "integrated.cca",
                           dims = 1:dims, reduction.name = reduction.name, 
                           verbose = TRUE, return.model = TRUE)
  print("<<<UMAP coordinates determined>>>")
  return(seurat.object)
}

generate_QC_plots <- function(seurat.object, reduction.name, target_expID, path) {
  # Ensure the output directory exists
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  p <- VlnPlot(seurat.object, features = c("nFeature_RNA", "nCount_RNA", 
                                           "percent.mt"), 
               ncol = 3, group.by = "expID", pt.size = 0)
  print(p)
  ggsave(filename = paste0(path, "/var_ftr_vln_expID.png"), plot = p)
  
  p <- VlnPlot(seurat.object, features = c("nFeature_RNA", "nCount_RNA", 
                                           "percent.mt"), ncol = 3, pt.size = 0)
  print(p)
  ggsave(filename = paste0(path, "/var_ftr_vln_per_cluster.png"), plot = p)
  
  p1 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA",
                       feature2 = "percent.mt", group.by = "expID")
  p2 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA",
                       feature2 = "nFeature_RNA", group.by = "expID")
  p <- p1 + p2
  print(p)
  ggsave(filename = paste0(path, "/var_ftr_scatter.png"), plot = p)
  
  #Plot regression variables on feature scatter plot
  ftrs = c("nCount_RNA", "percent.mt","percent.Hsp","nFeature_RNA")
  p <- FeaturePlot(seurat.object, features = ftrs, reduction = reduction.name, 
                   pt.size = 0.3)
  print(p)
  ggsave(filename = paste0(path, "/ftr_plot_regression_vars.png"), plot = p)
  
  #Plot the UMAPs split by experimental ID
  p <- DimPlot(seurat.object, reduction = reduction.name, label = FALSE,
               pt.size = 0.3, split.by = "expID")
  print(p)
  ggsave(filename = paste0(path, "/UMAP_unnamed_exp_split.png"), plot = p)
  
  #Plot a single target experiment vs the rest
  if(target_expID != ""){  
    target_subset <- subset(seurat.object, expID == target_expID)
    others_subset <- subset(seurat.object, expID != target_expID)
    p_target <- DimPlot(target_subset, reduction = reduction.name, 
                        label = FALSE, pt.size = 0.3) +
      ggtitle(paste("Target Experiment:", target_expID))
    p_other <- DimPlot(others_subset, reduction = reduction.name,
                       label = FALSE, pt.size = 0.3) +
      ggtitle("Other experiments, grouped")
    p <- p_target + p_other
    print(p)
    ggsave(filename = paste0(path, "/UMAP_unnamed_target_exp_v_rest.png"), plot = p)
  } else {
    message("target_expID is empty, skipping the plot.")
  }
  #Plot UMAP with overlayed experiments to check for batch effects
  colors <- brewer.pal(length(unique(seurat.object$expID)), "Set2") #Set colors
  p <- DimPlot(seurat.object, reduction = reduction.name, pt.size = 0.5,
               group.by = "expID", 
               cols = colors)
  print(p)
  ggsave(filename = paste0(path, "/UMAP_unnamed_exp_overlap.png"), plot = p)
} 

generate_ftr_plots_umap <- function(seurat.object, reduction.name, path) {
  # Ensure the output directory exists
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  #####
  #Feature plots & UMAP
  #####
  #Plot cell type discriminating genes on feature plot
  ftrs = c("CD14", "CD3E", "CD19", "NKG7", "CD4",
           "TNF", "IL6", "IL1B", "CD8A", "GZMA", "IL2RA", "NCAM1")
  p <- FeaturePlot(seurat.object, features = ftrs, reduction = reduction.name,
                   pt.size = 0.3)
  print(p)
  ggsave(filename = paste0(path, "/ftr_plots_genes.png"), plot = p, width = 16, height = 12)
  
  #####
  #UMAPs
  #####
  #Plot the combined UMAP
  p <- DimPlot(seurat.object, reduction = reduction.name, label = TRUE, 
               label.size = 5, combine = TRUE)
  print(p)
  ggsave(filename = paste0(path, "/UMAP_unnamed.png"), plot = p, width = 9, height = 5)
  
}

seurat_integrate_v2 <- function(
    #signature cell cycle
  signature_cell_cycle_human_mouse,
  #species
  species,
  cond.2.rmv, #Conditions to remove
  
  #UMI matrix
  matrix,
  #meta data
  meta_data,
  #Variables
  integrate_over, # Batch over "CC" or "expID"
  
  #Seurat
  min_cells_percent,
  min_gene_number,
  mito_cutoff,
  dims_use,
  #common variables to regress are: "n_UMI", "percent.mt", "Phase"
  vars_to_regress,
  join_layers,
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
    signature_cell_cycle <- subset(signature_cell_cycle_human_mouse, species == "human")
    mito_genes = "^MT-"
    ribo_genes = "^RPL"
    heatshock_genes = "^HSP"
  }
  if(species == "mouse"){
    signature_cell_cycle <- subset(signature_cell_cycle_human_mouse, species == "mouse")
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
  if(!is.null(cond.2.rmv)){
    seurat.object_v5 <- subset(seurat.object_v5, 
                               cells = which(!(seurat.object_v5@meta.data$expID %in% cond.2.rmv)))
    print("Experiments kept: ")
    print(unique(seurat.object_v5@meta.data$expID))
  }
  
  #Extract factor over which to integrate
  if(!is.null(integrate_over)){
    integrate_by <- seurat.object_v5[[integrate_over]]
    integrate_by_FCI <- rownames(integrate_by)
    integrate_by_type <- integrate_by[,1]
    names(integrate_by_type) <- integrate_by_FCI
    
    #Split data
    seurat.object_v5[["RNA"]] <- split(seurat.object_v5[["RNA"]], 
                                       f = integrate_by_type)
    print("<<<Seurat layers according to integration parameter set>>>")
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
  DefaultLayer(seurat.object_v5[["RNA"]]) <- 'data'
  
  #Find variable genes, scale and PCA
  seurat.object_v5 <- FindVariableFeatures(seurat.object_v5, verbose = FALSE)
  
  seurat.object_v5 <- ScaleData(seurat.object_v5,
                                assay = "RNA",
                                model.use = "linear",
                                features = rownames(seurat.object_v5@assays$RNA),
                                vars.to.regress = vars_to_regress,
                                verbose = TRUE)
  
  seurat.object_v5 <- RunPCA(seurat.object_v5, npcs = dims_use, verbose = FALSE)
  
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
  
  stdevs <- seurat.object_v5@reductions$pca@stdev
  min_stdev <- min(stdevs)
  thresh <- min_stdev * 1.05 # Make a st_dev threshold @ 5% deviation
  #Return the dimension which fulfills the above threshold criterium, assign value to parameters
  dims_use <- which(stdevs <= thresh)[1]
  cat("PCA dimensions set: ", dims_use, "\n")
  print("<<<FVB, Scale, Regression, PCA done>>>")
  
  #Compute Integrations
  # Ensure k.weight is smaller than the minimum number of cells
  k.weight <- floor(min(100, min_cells))
  cat("k.weight set to ", k.weight, "\n")
  
  #Check number of count matrices in Layers
  layer_names <- names(seurat.object_v5@assays$RNA@layers)
  count_layers <- grep("count", layer_names, value = TRUE)
  
  if(length(count_layers) > 1){
    seurat.object_v5 <- IntegrateLayers(
      object = seurat.object_v5, method = CCAIntegration,
      orig.reduction = "pca", new.reduction = "integrated.cca",
      k.weight = k.weight,
      verbose = FALSE)
    print("<<<Multiple layers integrated>>>")
  }
  
  #Rename the "pca" reduction as "integrated.cca" if there is only 1 layer
  if(length(seurat.object_v5@reductions) == 1 && "pca" %in% names(seurat.object_v5@reductions)){
    # Rename the reduction
    seurat.object_v5@reductions[["integrated.cca"]] <- seurat.object_v5@reductions[["pca"]]
    seurat.object_v5@reductions[["pca"]] <- NULL
    
    print("Reduction renamed from 'pca' to 'integrated.cca'.")
  } else {
    print("Renaming not performed. Either there are multiple reductions or 'pca' does not exist.")
  }
  
  if(join_layers == "yes"){
    seurat.object_v5 <- JoinLayers(seurat.object_v5)
    print("<<<Layers joined>>>")
  }
  
  #Cell cycle scoring
  s.genes <- subset(signature_cell_cycle, cycle_state == "S")$GeneSymbol
  g2m.genes <- subset(signature_cell_cycle, cycle_state == "G2/M")$GeneSymbol
  
  seurat.object_v5 <- CellCycleScoring(seurat.object_v5,
                                       s.features = s.genes,
                                       g2m.features = g2m.genes)
  print("<<<Cell cycle score computed>>>")
  
  #Find neighbors
  seurat.object_v5 <- FindNeighbors(seurat.object_v5, reduction = "integrated.cca", 
                                    dims = 1:dims_use, verbose = TRUE)
  print("<<<Neighbors found>>>")
  return(list(seurat.object_v5, dims_use))
}