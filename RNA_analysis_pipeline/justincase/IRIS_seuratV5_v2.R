#Author: Joern Pezoldt & T Ferrari
#Date init: 12.01.2024
#Scopes:
## 2) Superimposte image data classification for DEG analysis and visualization

"Note: all combined integrations are CCA integrated by layer"
# install.packages("languageserver")
options(future.globals.maxSize = +Inf)
#libraries
library(dplyr)
library(Seurat)
library(SeuratData)
library(UpSetR)

library(topGO)
library(fgsea)
library(reticulate)
library(leidenAlg)
library(openxlsx)
library(harmony)

library(grid)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(scales)
library(ggrepel)
library(EnhancedVolcano)
library(stringr)

#load custom functions around scRNA-seq with IRIS, imaging and tools used for analysis
PATH_functions <- "/home/tferrari/updepla/users/tferrari/GitHub_repos/IRIS_Analysis_pipelines/RNA_analysis_pipeline/PBMC_seurat_functions"
source(paste0(PATH_functions,"/seuratV5_integration_v2.R"))
source(paste0(PATH_functions,"/support_transcriptome_integration.R"))
source(paste0(PATH_functions, "/count_annotation_matrix_creation.R"))
source(paste0(PATH_functions, "/create_and_save_plots.R"))
source(paste0(PATH_functions, "/calculate_signatures_markers.R"))
source(paste0(PATH_functions, "/rename_clusters_calculate_GO.R"))
source(paste0(PATH_functions, "/create_annotation_dataframe.R"))

#####
#Set PATHs
#####
#UserID
userID <- "tferrari"
analysisID <- "2024_PBMCs"
#Samples to integrate
sample_ids_IRIS <-   c("JP303", "TF050", "JP304", # KMT2Ar no blasts (healthy profile)
                       "TF047", "TF048", "NG079", "JP318", "TF059", # KMT2Ar PB samples
                       "TF062", "JP316")# , # "NG084", # KMT2Ar BM
  # c("NG037", "NG039", "NG042", "NG050", "NG078", "NG075", "JP314", "TF061", "TF057") # HD PB controls

  # c("NG064", "TF051", "TF052", "JP317", "TF060", "TF064", "TF066", "JP315", "TF063", "NG083", # NPM1 PB
  #                      "TF058", "NG080", "NG082", "TF065") # NPM1 BM
sample_ids <- c(sample_ids_IRIS)
#Define paths
PATH_input_IRIS_sequencing <- paste0("/home/tferrari/updepla/projects/iris/4_sequencing_analysis")
PATH_input_IRIS_imaging <- paste0("/home/tferrari/updepla/users/tferrari/Experiments/Image_Analysis/AML_experiments/experiment_files")

PATH_output <- paste0("/home/",userID,
                      "/updepla/users/", userID, 
                      "/Experiments/RNA_analysis/IRIS/AML_experiments/results/AML_samples/KMT2Ar_BM-PB/harmony")
PATH_output_figures <- paste0(PATH_output,"/plots")
PATH_output_objects <- paste0(PATH_output,"/objects")
PATH_output_tables <- paste0(PATH_output,"/tables")
l.result.paths <- list(PATH_output_tables, PATH_output_objects, PATH_output_figures)

for(path in l.result.paths){
  if(!file.exists(path)){
    dir.create(path, recursive = TRUE)
    message(paste("Directory created: ", path))
  } else {
    message(paste("Directory already exists: ", path))
  }
}

########
#Signatures
########
#Path signatures
PATH_signature <- "/home/tferrari/updepla/users/tferrari/Experiments/RNA_Analysis/IRIS/AML_experiments/signatures"
#Signatures PBMCs
signatures.PBMCs <- read.delim(paste0(PATH_signature,"/healthy_sig_top50.txt"),
                               sep = "\t", header = TRUE)
signatures.PBMCs.250degs <- read.delim(paste0(PATH_signature,"/healthy_sig_top250.txt"),
                                       sep = "\t", header = TRUE)
signatures.PBMCs.500degs <- read.delim(paste0(PATH_signature,"/healthy_sig_top500.txt"),
                                       sep = "\t", header = TRUE)
signatures.BMMCs <- read.delim(paste0(PATH_signature,"/NPM1_NG064_BMMC_TSP21_TSP25_sig.txt"),
                               sep = "\t", header = TRUE)
# Signatures cell cycle
signature_cell_cycle_human_mouse <- readRDS(paste0(PATH_signature,
                                                   "/signature_cell_cycle_human_mouse.rds"))

#####
#Parameters to be modified
#####
#Variables
species = "human"
# Batch over "CC" or "expID"
integrate_over = "expID" # Next should batch by "sampleID" and "cohort"

#Threshing
min_cells_percent = 0.005
min_gene_number = 500
mito_cutoff = 20
min_nUMI = 2500

#Set starting dimensions
dims_use = 75
#common variables to regress are: "n_UMI", "percent.mt", "Phase"
vars_to_regress <- c()

set.seed(300)

#Set threshold values for marker DEGs
min.pct <- 0.1
logfc.threshold <- 0.25

#Set experiment IDs to be removed if necessary (only necessary if 
#want to maintain "uniqueness" of clusters for datasets with low relative cell counts)
cond.2.rmv <- c()

#Set UMAP reduction name
#Set whether to first split layers and then perform IntegrateLayers on the seurat.object
#Will split and integrate layers if "yes"
split.layers <- "yes"
integrate.layers <- "yes"
integration <- "pca"
reduction <- "harmony" # Option: "cca", "rpca", "harmony"
#Select the type of community clustering algorithm
cluster_algo <- "Leiden"

#######################################
#Image data
#######################################
#####
# Load Count matrix IRIS platform
#####
#Annotate doublets and store them in a list
l.image.doublet.ROI.annotation <- annotate_doublets(PATH_input_IRIS_imaging, 
                                                    sample_ids_IRIS)
#LiquideDrop data
l.meta.mtx <- create_matrix_meta_data(sample_ids_IRIS, 
                                      l.image.doublet.ROI.annotation)
#####
#Batch integration with Seurat
#####
#Seurat
print(unique(l.meta.mtx[[1]]$CellType_Condition))
print(unique(l.meta.mtx[[1]]$expID))

integration.results <- seurat_integrate(
  signature_cell_cycle_human_mouse,
  species,
  cond.2.rmv,
  l.meta.mtx[[2]],
  l.meta.mtx[[1]],
  integrate_over,
  min_cells_percent,
  min_gene_number,
  mito_cutoff,
  dims_use,
  vars_to_regress,
  split.layers,
  integrate.layers,
  reduction,
  paste0(PATH_output_figures, "/QC_plots"))

seurat.object <- integration.results[[1]]

# saveRDS(file = paste0(PATH_output_objects, "/seurat_object_no_neighbors.rds"),
        # seurat.object)
# seurat.object <- readRDS(file = paste0(PATH_output_objects, "/seurat_object_no_neighbors.rds"))

dims_use <- 58
#Find neighbors
seurat.object <- FindNeighbors(seurat.object, reduction = reduction, 
                               dims = 1:dims_use, verbose = TRUE)
print("<<<Neighbors found>>>")

saveRDS(file = paste0(PATH_output_objects, "/seurat_object_neighbors.rds"),
        seurat.object)

# seurat.object <- readRDS(file = paste0(PATH_output_objects, "/seurat_object_neighbors.rds"))

#########################
# Plot dimensional and resolution QC plots
#########################
#######
# Plot "low resolution" cluster for global visualization
#######
cluster_algo <- "Leiden"

print(paste("Algorithm used: ", cluster_algo))
if(cluster_algo == "Leiden"){
  algo = 4
} else if(cluster_algo == "Louvain") {
  algo = 1
}

#Plot the UMAPs split by experimental ID once so that we can appreciate expID splits, res unimportant
seurat.object_qc <- FindClusters(seurat.object, resolution = 0.7, algorithm = algo,
                                 verbose = TRUE)
#Run UMAP with 0.7 "standard" resolution, just for overview
seurat.object_qc <- RunUMAP(seurat.object_qc, dims = 1:dims_use, 
                            reduction = reduction,  verbose = TRUE)

p <- DimPlot(seurat.object_qc, reduction = "umap", label = FALSE,
             pt.size = 0.5, split.by = "expID", ncol = 4) + NoLegend()
print(p)
ggsave(
  filename = paste0(PATH_output_figures, "/QC_plots", "/UMAP_exp_split.pdf"), 
  plot = p,
  width = 16, 
  height = ceiling(length(unique(seurat.object_qc@meta.data$expID)) / 4) * 4, 
  limitsize = FALSE, device = 'pdf', units = "in", useDingbats = FALSE
)

#Plot UMAP with overlayed experiments to check for batch effects
# colors <- brewer.pal(length(unique(seurat.object$expID)), "Set2") #Set colors
p <- DimPlot(seurat.object_qc, reduction = "umap", pt.size = 0.5,
             group.by = "expID") + #, cols = cluster_colors) +
  print(p)
ggsave(filename = paste0(PATH_output_figures, "/QC_plots", "/UMAP_exp_overlap.pdf"), plot = p,
       width = 16, height = 8, device = 'pdf', units = "in", useDingbats = FALSE)

########
# Signature z-score plotting
########
# Temprorarily plot and calculate signature z-scores to decide on resolutions
#Select top 500 DEGs per cluster and order in descending for avg_log2FC and plot
seurat.object_qc <- JoinLayers(seurat.object_qc)

"All quoted below is necessary for the calculation of AddModuleScore, the precision of 
which is not optimal"
# Select appropriate number of DEGs for signature printing
signatures.PBMCs %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  dplyr::filter(avg_log2FC > 0.25) %>%
  slice_head(n = 100) %>%
  ungroup() -> top.PBMC.DEGs
# Enforce underscore as spacer
top.PBMC.DEGs$cluster <- gsub(" ", "_", top.PBMC.DEGs$cluster)
top.PBMC.DEGs$cluster <- paste0(top.PBMC.DEGs$cluster, "_(PBMC)")

signatures.BMMCs %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  dplyr::filter(avg_log2FC > 0.25) %>%
  slice_head(n = 100) %>%
  ungroup() -> top.BMMC.DEGs
# Add suffix for overlapping signatures
top.BMMC.DEGs$cluster <- paste0(top.BMMC.DEGs$cluster, "_(BMMC)")

gene_list <- list()
for(score_name in unique(top.PBMC.DEGs$cluster)){
  print(score_name)
  gene_list[[score_name]] <- top.PBMC.DEGs[top.PBMC.DEGs$cluster == (score_name), ]$gene
  seurat.object_qc <- AddModuleScore(seurat.object_qc,
                                     features = list(gene_list[[score_name]]),
                                     name = score_name)
}
for(score_name in unique(top.BMMC.DEGs$cluster)){
  print(score_name)
  gene_list[[score_name]] <- top.BMMC.DEGs[top.BMMC.DEGs$cluster == (score_name), ]$gene
  seurat.object_qc <- AddModuleScore(seurat.object_qc, features = list(gene_list[[score_name]]),
                                     name = score_name)
}

# Remove the "1" suffix from the new score columns in meta.data
# Renaming columns in the meta.data slot
names(seurat.object_qc@meta.data)[(ncol(seurat.object_qc@meta.data) - length(gene_list) + 1):ncol(seurat.object_qc@meta.data)] <- names(gene_list)

# Remove underscores from gene names for plotting
features_clean <- gsub("_", " ", names(gene_list))  # or replace "_" with " " for spaces

mincoff <- c("q60")
plots <- FeaturePlot(seurat.object_qc, features = names(gene_list),
                     reduction = "umap", cols = c("lightgrey", "deepskyblue", "blue"),
                     min.cutoff = mincoff,
                     max.cutoff = c("q99"),
                     pt.size = 0.8, order = TRUE) 
# Update the titles of the individual plots to the cleaned gene names
plots <- lapply(seq_along(plots), function(i) {
  plots[[i]] + ggplot2::labs(title = features_clean[i]) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5,
                                margin = margin(b = 15)),
      plot.margin = margin(20, 20, 30, 20)
    )
})
p <- wrap_plots(plots) + plot_annotation(title = "Feature z-scores") &
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) & plot_layout(guides = 'collect')
print(p)
ggsave(filename = paste0(PATH_output_figures, "/z-score_plots",
                         "/pbmc_top50_z-score_mincoff_", mincoff, ".pdf"),
       plot = p, width = 32, height = 32, device = 'pdf', units = "in", useDingbats = FALSE)

if(!file.exists(paste0(PATH_output_figures, "/QC_plots/UMAP_plots"))){
  dir.create(paste0(PATH_output_figures, "/QC_plots/UMAP_plots"), recursive = TRUE)
  message(paste("Directory created: ", paste0(PATH_output_figures, "/QC_plots/UMAP_plots")))
} else {
  message(paste("Directory already exists: ", paste0(PATH_output_figures, "/QC_plots/UMAP_plots")))
}

#######
# UMAP resolution plotting for adequate cluster discrimination
######
resolutions <- seq(0.5, 3.0, by = 0.1) # Define the range of resolutions
for(res in resolutions){
  # Run clustering with different resolutions
  seurat.object_umapqc <- FindClusters(seurat.object, resolution = res, algorithm = algo,
                                       verbose = TRUE)
  #Run UMAP with current resolution
  seurat.object_umapqc <- RunUMAP(seurat.object_umapqc, dims = 1:dims_use, 
                                  reduction = reduction,
                                  verbose = TRUE)
  plot <- DimPlot(seurat.object_umapqc, label = TRUE, pt.size = 0.3) +
    labs(title = paste("UMAP clustering at resolution", res)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  print(plot)
  ggsave(filename = paste0(PATH_output_figures, "/QC_plots/UMAP_plots", "/umap_res", gsub("\\.", "_", res), ".pdf"), 
         plot = plot, height = 8, width = 10, device = 'pdf', units = "in", useDingbats = FALSE)
}
# 
#   seurat.object_v5 <- JoinLayers(seurat.object_v5)
#   markers.all <- FindAllMarkers(seurat.object_v5, only.pos = FALSE,
#                                 min.pct = min.pct, logfc.threshold = logfc.threshold)
# 
#   markers.all <- markers.all %>% arrange(cluster, desc(avg_log2FC))
#   head(markers.all)
# 
#   #Select top 250 DEGs
#   markers.all %>% group_by(cluster) %>%
#     top_n(n = 7, wt = avg_log2FC) -> top7
# 
#   #Plot heatmap
#   ordered_palette <- custom.palette[order(as.integer(names(custom.palette)))]
#   p <- DoHeatmap(seurat.object_v5, features = top7$gene,
#                  size = 3, angle = 90) +
#     # scale_fill_gradient2( low = rev(c('#0000CD','#000080', '#00008B')),
#     #                       mid = "white", high = rev(c('#FFD700','#EEC900','#CDAD00')),
#     #                       midpoint = 0, guide = "colourbar", aesthetics = "fill") +
#     NoLegend()
#   print(p)
#   ggsave(filename = paste0(path, "/deg_res", gsub("\\.", "_", res), ".pdf"),
#          plot = p, width = 10, height = 18, device = 'pdf', units = "in", useDingbats = FALSE)
# }

if(!file.exists(paste0(PATH_output_figures, "/unnamed_plots"))){
  dir.create(paste0(PATH_output_figures, "/unnamed_plots"), recursive = TRUE)
  message(paste("Directory created: ", paste0(PATH_output_figures, "/unnamed_plots")))
} else {
  message(paste("Directory already exists: ", paste0(PATH_output_figures, "/unnamed_plots")))
}

res <- 2.2

seurat.object <- cluster_UMAP_seurat(seurat.object, res, dims_use,
                                     reduction, cluster_algo, paste0(PATH_output_figures, "/unnamed_plots"))

# Prepare the text to be printed
text_content <- paste(
  "##### Parameters set for refinement and integration #####\n",
  "sample_IDs: ", paste(sample_ids_IRIS, collapse = ", "), "\n",
  "Species: ", species, "\n",
  "Integrate Over: ", integrate_over, "\n",
  "Minimum Cells Percent: ", min_cells_percent, "\n",
  "Minimum Gene Number: ", min_gene_number, "\n",
  "Mitochondrial Cutoff: ", mito_cutoff, "%\n",
  "Minimum nUMI: ", min_nUMI, "\n",
  "Dimensions to Use: ", dims_use, "\n",
  "Variables to Regress: ", paste(vars_to_regress, collapse = ", "), "\n",
  "Seed Value: ", 300, "\n",
  "Minimum Percentage for Marker DEGs: ", min.pct, "\n",
  "Log Fold Change Threshold: ", logfc.threshold, "\n",
  "Conditions to Remove: ", paste(cond.2.rmv, collapse = ", "), "\n",
  "Split Layers: ", split.layers, "\n",
  "Integrate Layers: ", integrate.layers, "\n",
  "Clustering Algorithm: ", cluster_algo, "\n",
  "UMAP Resolution: ", res, "\n",
  "UMAP Reduction Name: ", reduction, "\n",
  "Notes: "
)

writeLines(text_content, con = paste0(PATH_output, "/integration_parameters.txt"))

#####
# Signature plotting signatures
#####
#Select top 500 DEGs per cluster and order in descending for avg_log2FC and plot
seurat.object <- JoinLayers(seurat.object) # Also necessary for DEG calculation

#Save seurat.object
saveRDS(file = paste0(PATH_output_objects, "/seurat_object_integrated.rds"),
        seurat.object)

print(paste("Fully integrated seurat object saved to: ",
            paste0(PATH_output_objects, "/seurat_object_integrated.rds")))

# seurat.object <- readRDS(file = paste0(PATH_output_objects, "/seurat_object_integrated.rds"))

# "All quoted below is necessary for the calculation of AddModuleScore, the precision of 
# which is not optimal"
# # Select appropriate number of DEGs for signature printing
# signatures.PBMCs %>%
#   group_by(cluster) %>%
#   arrange(desc(avg_log2FC), .by_group = TRUE) %>%
#   dplyr::filter(avg_log2FC > 0.25) %>%
#   slice_head(n = 100) %>%
#   ungroup() -> top.PBMC.DEGs
# # Enforce underscore as spacer
# top.PBMC.DEGs$cluster <- gsub(" ", "_", top.PBMC.DEGs$cluster)
# top.PBMC.DEGs$cluster <- paste0(top.PBMC.DEGs$cluster, "_(PBMC)")
# 
# signatures.BMMCs %>%
#   group_by(cluster) %>%
#   arrange(desc(avg_log2FC), .by_group = TRUE) %>%
#   dplyr::filter(avg_log2FC > 0.25) %>%
#   slice_head(n = 100) %>%
#   ungroup() -> top.BMMC.DEGs
# # Add suffix for overlapping signatures
# top.BMMC.DEGs$cluster <- paste0(top.BMMC.DEGs$cluster, "_(BMMC)")
# 
# gene_list <- list()
# for(score_name in unique(top.PBMC.DEGs$cluster)){
#   print(score_name)
#   gene_list[[score_name]] <- top.PBMC.DEGs[top.PBMC.DEGs$cluster == (score_name), ]$gene
#   seurat.object <- AddModuleScore(seurat.object,
#                                   features = list(gene_list[[score_name]]),
#                                   name = score_name)
# }
# for(score_name in unique(top.BMMC.DEGs$cluster)){
#   print(score_name)
#   gene_list[[score_name]] <- top.BMMC.DEGs[top.BMMC.DEGs$cluster == (score_name), ]$gene
#   seurat.object <- AddModuleScore(seurat.object, features = list(gene_list[[score_name]]),
#                                   name = score_name)
# }
# 
# # Remove the "1" suffix from the new score columns in meta.data
# # Renaming columns in the meta.data slot
# names(seurat.object@meta.data)[(ncol(seurat.object@meta.data) - length(gene_list) + 1):ncol(seurat.object@meta.data)] <- names(gene_list)
# 
# # Remove underscores from gene names for plotting
# features_clean <- gsub("_", " ", names(gene_list))  # or replace "_" with " " for spaces
# 
# mincoff <- c("q60")
# plots <- FeaturePlot(seurat.object, features = names(gene_list),
#                      reduction = "umap", cols = c("lightgrey", "deepskyblue", "blue"),
#                      min.cutoff = mincoff,
#                      max.cutoff = c("q99"),
#                      pt.size = 0.8, order = TRUE) 
# # Update the titles of the individual plots to the cleaned gene names
# plots <- lapply(seq_along(plots), function(i) {
#   plots[[i]] + ggplot2::labs(title = features_clean[i]) +
#     theme(
#       plot.title = element_text(size = 12, face = "bold", hjust = 0.5,
#                                 margin = margin(b = 15)),
#       plot.margin = margin(20, 20, 30, 20)
#     )
# })
# p <- wrap_plots(plots) + plot_annotation(title = "Feature z-scores") &
#   theme(plot.title = element_text(hjust = 0.5, face = "bold")) & plot_layout(guides = 'collect')
# print(p)
# ggsave(filename = paste0(PATH_output_figures, "/unnamed_plots",
#                          "/ct_ftr_mincoff_", mincoff, ".pdf"),
#        plot = p, width = 32, height = 32, device = 'pdf', units = "in", useDingbats = FALSE)
# 
# saveRDS(file = paste0(PATH_output_objects, "/seurat_object_zscore.Rds"), seurat.object)
# print(paste("seurat object with Z-scores in @meta_data saved to: ",
#             paste0(PATH_output_objects, "/seurat_object_zscore.Rds")))

# seurat.object <- readRDS(file = paste0(PATH_output_objects, "/seurat_object_zscore.Rds"))

#####
#DEGs general
#####
l.markers <- calculate_and_plot_marker_DEGs(seurat.object, min.pct,
                                            logfc.threshold, 1000,
                                            paste0(PATH_output_figures,
                                                   "/unnamed_plots/DEG_heatmap.pdf"))
# Also save it to an excel file for easy reading
wb <- createWorkbook() #Create new workbook
addWorksheet(wb, "Unnamed_all_markers") #Add worksheet
writeData(wb, "Unnamed_all_markers", l.markers[[1]]) #Write data
addWorksheet(wb, "Unnamed_top_markers") #Add worksheet
writeData(wb, "Unnamed_top_markers", l.markers[[2]]) #Write data
addWorksheet(wb, "Unnamed_bottom_markers") #Add worksheet
writeData(wb, "Unnamed_bottom_markers", l.markers[[3]]) #Write data
saveWorkbook(wb, paste0(PATH_output_tables, "/unnamed_markers.xlsx")) #Save workbook

#############################
#Plotting
#############################
#####
#QC
#####
generate_QC_plots(seurat.object, res, "umap", c(head(sample_ids_IRIS, -1)), 
                  paste0(PATH_output_figures, "/QC_plots"))
print(paste("QC plots saved to: ", paste0(PATH_output_figures, "/QC_plots")))

generate_bar_plots(seurat.object, paste0(PATH_output_figures, "/QC_plots"))
print(paste("Barplots saved to: ", paste0(PATH_output_figures, "/QC_plots")))

#####
# Plot the UMAP with binned cohort colors
#####
kmt2ar.pb <- c("TF047", "TF048", "NG079", "JP318", "TF059")  # KMT2Ar PB samples
kmt2ar.bm <- c("TF062", "JP316")
kmt2ar.hd.pb <- c("JP303", "TF050", "JP304")

npm1.pb <- c("NG064", "TF051", "TF052", "JP317", "TF064", "TF066", "JP315", "TF063", "NG083")
npm1.bm <- c("TF058", "NG080", "NG082", "TF065")

hd.pb <- c("NG037", "NG039", "NG042", "NG050", "NG078", "NG075", "JP314", "TF061", "TF057")

all_expIDs <- unique(seurat.object@meta.data$expID)

expID.colors <- rep(alpha("lightgray", 0.1), length(all_expIDs))
names(expID.colors) <- all_expIDs

expID.colors[kmt2ar.pb] <- rep(alpha("orange", 1.0))
expID.colors[kmt2ar.bm] <- rep(alpha("red"), 1.0)
expID.colors[kmt2ar.hd.pb] <- rep(alpha("green"), 1.0)
expID.colors[npm1.pb] <- rep(alpha("cyan"), 1.0)
expID.colors[npm1.bm] <- rep(alpha("darkblue"), 1.0)

if(!file.exists(paste0(PATH_output_figures, "/QC_plots/cohort_plots"))){
  dir.create(paste0(PATH_output_figures, "/QC_plots/cohort_plots"), recursive = TRUE)
  message(paste("Directory created: ", paste0(PATH_output_figures, "/QC_plots/cohort_plots")))
} else {
  message(paste("Directory already exists: ", paste0(PATH_output_figures, "/cohort_plots")))
}
p <- DimPlot(seurat.object, group.by = "expID", cols = expID.colors)
print(p)
ggsave(filename = paste0(PATH_output_figures, "/QC_plots/cohort_plots/UMAP_cohort_overlap.pdf"),
       height = 12, width = 14, device = 'pdf', units = "in", useDingbats = FALSE)

all_cohorts <- list(kmt2ar.pb, kmt2ar.bm, kmt2ar.hd.pb)
  list(hd.pb)
  list(npm1.pb, npm1.bm) 
#  hd.pb) #, bm.osc)
cohort_names <- c("kmt2ar_pb", "kmt2ar_bm", "kmt2ar_hd_pb")
  c("hd.pb")
  c("npm1_pb", "npm1_bm")
                  # "hd_pb") #, "bm_osc")

i = 1
for(cohort in all_cohorts){
  print(cohort)
  sub.seurat <- subset(seurat.object, 
                       subset = seurat.object@meta.data$expID %in% cohort)
  Idents(sub.seurat) <- sub.seurat@meta.data$seurat_clusters
  plot_df <- FetchData(sub.seurat, 
                       vars = c("umap_1", "umap_2", "expID", "seurat_clusters"))
  
  centroids <- aggregate(cbind(umap_1, umap_2) ~ seurat_clusters, 
                         data = plot_df, FUN = mean)
  
  p <- ggplot(plot_df, aes(x = umap_1, y = umap_2, color = expID)) +
    geom_point(size = 0.5, alpha = 1.0) +
    scale_color_manual(values = expID.colors) +
    geom_text_repel(data = centroids,
                    aes(x = umap_1, y = umap_2, label = seurat_clusters),
                    color = "black", size = 4, fontface = "bold",
                    box.padding = 0.4, point.padding = 0.4) +
    # Set global axis limits
    # coord_cartesian(xlim = global_xlim, ylim = global_ylim) +
    theme_minimal() +
    labs(color = "Experiment ID") +
    theme(
      legend.position = "right"
    )
  print(p)
  ggsave(filename = paste0(PATH_output_figures, "/QC_plots/cohort_plots/", cohort_names[[i]], "_UMAP.pdf"),
         height = 8, width = 9, device = 'pdf', units = "in", useDingbats = FALSE)
  i=i+1
}

########
# Feature plots
########
plot_specific_ftr(seurat.object, paste0(PATH_output_figures, "/feature_plots"))
print(paste("Class specific feature plots saved to ",
            paste0(PATH_output_figures, "/feature_plots")))

########
# Plots per sample ID
########
# First create sample_id variable in meta.data and then assign values for updepla_code (is a sample ID but sampleID is already taken)
# depending on expID value, also provide location origin afterwards (will be useful in multiple group.by)
seurat.object@meta.data <- seurat.object@meta.data %>% 
  mutate(updepla_code = case_when(
    expID %in% c("JP303", "JP304", "TF050") ~ "a_CH5_AML",
    expID == "TF047" ~ "a_CH6_AML",
    expID == "TF048" ~ "a_CH8_AML",
    expID %in% c("JP318", "TF062") ~ "a_SL1_AML",
    expID %in% c("TF059", "NG079", "JP316") ~ "a_SL2_AML",
    TRUE ~ NA_character_
  ))

seurat.object@meta.data <- seurat.object@meta.data %>% 
  mutate(location = case_when(
    expID %in% c("JP303", "JP304", "TF050") ~ "peripheral_blood",
    expID == "TF047" ~ "peripheral_blood",
    expID == "TF048" ~ "peripheral_blood",
    expID == "JP318" ~ "peripheral_blood",
    expID == "TF062" ~ "bone_marrow",
    expID %in% c("TF059", "NG079") ~ "peripheral_blood",
    expID == "JP316" ~ "bone_marrow",
    TRUE ~ NA_character_
  ))

seurat.object@meta.data <- seurat.object@meta.data %>% 
  mutate(primary_mut = case_when(
    expID %in% c("JP303", "JP304", "TF050") ~ "t(10;11)(p12;q23.3);KMT2A::MLLT10",
    expID == "TF047" ~ "t(11;19)(q23.3;p13.3);KMT2A::MLLT1",
    expID == "TF048" ~ "t(6;11)(q27;q23.3);KMT2A::AFDN",
    expID %in% c("JP318", "TF062") ~ "t(11;19)(q23;p13);KMT2A::ELL",
    expID %in% c("TF059", "NG079", "JP316") ~ "t(10;11)(p12;q23);KMT2A::MLLT10",
    TRUE ~ NA_character_
  ))

# Combine mutation and location for plotting
seurat.object@meta.data$pr_mut_loc <- paste(seurat.object$primary_mut, seurat.object$location)

saveRDS(file = paste0(PATH_output_objects, "/seurat_object_metadata.rds"), seurat.object)
print(paste("Final fully integrated seurat object saved to: ",
            paste0(PATH_output_objects, "/seurat_object_metadata.rds")))

seurat.object <- readRDS(paste0(PATH_output_objects, "/seurat_object_metadata.rds"))

p <- DimPlot(seurat.object, split.by = c("updepla_code"), ncol = 3) +
  theme(strip.text = element_text(size = 8))
print(p)
ggsave(filename = paste0(PATH_output_figures, '/unnamed_plots/split_umaps/umap_updepla_code.pdf'),
       plot = p, width = 14, 
       height = ceiling(length(unique(seurat.object@meta.data$updepla_code))/3)*4,
       device = 'pdf', units = "in", useDingbats = FALSE)

######
#Cluster annotation
######
#Rename clusters with major cell group attribution and then subselect again
#Do not use the same names as for signature_PBMCs
new.cluster.names <- c("", "", # Cluster 1-2
                       "", "CD14+_classical_monocytes", # Cluster 3-4
                       "HSPC-like_myeloid_blast_(CD34h_CD117h_", "", # Cluster 5-6
                       "", "", # Cluster 7-8
                       "", "", # Cluster 9-10
                       "HSPC-like_myeloid_blast_(CD34h_CD117h_", "", # Cluster 11-12
                       "", "", # Cluster 13-14
                       "", "", # Cluster 15-16
                       "", "", # Cluster 17-18
                       "CD16+_non-classical_monocytes", "", # Cluster 19-20
                       "CD19+_B_cells", "", # Cluster 21-22
                       "HSPC-like_myeloid_blast_(CD34h_CD117h_", "" # Cluster 23-24
                       )

seurat.object <- rename_clusters(seurat.object, new.cluster.names,
                                 paste0(PATH_output_figures, "/named_plots"))

saveRDS(file = paste0(PATH_output_objects, "/seurat_object_final.rds"), seurat.object)
print(paste("Final fully integrated seurat object saved to: ",
            paste0(PATH_output_objects, "/seurat_object_final.rds")))

# seurat.object <- readRDS(file = paste0(PATH_output_objects, "/seurat_object_final.Rds"))

####
# Plot the DNA microarray, PCR, FACS clinical profile
###
# Subset seurat object if necessary
sample_path <- paste0(PATH_output_figures, "/sample_subset_plots", "/a_CH5_AML")
if(!file.exists()){
  dir.create(sample_path, recursive = TRUE)
  message(paste("Directory created: ", sample_path))
} else {
  message(paste("Directory already exists: ", sample_path))
}

p <- FeaturePlot(seurat.object, features = c(),
                 ncol = 4) #Non-classical monocytes
print(p)
ggsave(filename = paste0(sample_path, "/"), 
       plot = p, width = 16, height = 20, device = 'pdf', units = "in", useDingbats = FALSE)

##########
# Detect new markers
##########
l.new.markers <- calculate_and_plot_marker_DEGs(seurat.object, min.pct, logfc.threshold, 
                               1000, paste0(PATH_output_figures,
                                            "/named_plots/named_DEG_heatmap.pdf"))
# new.markers.all <- l.new.markers[[1]]
# topn <- l.new.markers[[2]]
# bottomn <- l.new.markers[[3]]

#Also save it to an excel file for easy reading
wb <- createWorkbook() #Create new workbook
addWorksheet(wb, "named_markers_all") #Add worksheet
writeData(wb, "named_markers_all", l.new.markers[[1]]) #Write data
addWorksheet(wb, "named_top_markers") #Add worksheet
writeData(wb, "named_top_markers", l.new.markers[[2]]) #Write data
addWorksheet(wb, "named_bottom_markers") #Add worksheet
writeData(wb, "named_bottom_markers", l.new.markers[[3]]) #Write data
saveWorkbook(wb, paste0(PATH_output_tables, "/named_markers.xlsx")) #Save workbook

saveRDS(l.new.markers, file = paste0(PATH_output_tables, "/l_markers.rds"))
print(paste("Markers list saved to: ", paste0(PATH_output_tables, "/l_markers.rds")))
print("Note: Markers list 1: All markers, 2: top 500 EDGs, 3: bottom 500 DEGs")

# l.new.markers <- readRDS(file = paste0(PATH_output_tables, "/l_markers.Rds"))

########
# Violin plots for specific blast/oncological markers
########

#Find genes with specific sequences
# Extract all gene names from the Seurat object
# all_genes <- rownames(seurat.object)
# target_genes <- grep("HB", all_genes, value = TRUE, ignore.case = TRUE)

mrkrs.2.plot <- list("KIT",
                     "CD8A", "GZMB", "PRF1", # "CD8+ Cytotoxic T cells"
                     "CD14", "FCGR3A", "LYZ", "S100A9", # FCGR3A (low/–), "CD14+ Classical Monocytes"
                     "CD4", "CCR7", "SELL", # "CD4+ Naive T cells"
                     "IL7R", # "CD4+ Central Memory T cells"
                     "PTPRC", "GZMK", # "CD8+ Effector Memory T cells"
                     "NCAM1", "NKG7", "GNLY", "KLRD1", # "CD56+ Natural Killer cells"
                     "CD19", "MS4A1", "CD79A", "CD79B", "CD22", # "CD19+ B cells"
                     "FOXP3", "IL2RA", "CTLA4", # "CD4+ Regulatory T (T-reg) cells"
                     "KLRB1", "SLC4A10", "IL18RAP",# "Mucosal-Associated Invariant T (MAIT) cells"
                     "GATA3", "CCR4", "IL4R", # "T Helper 2 (Th2) T-reg cells"
                     "S100A8", # CD14 (low/–); "CD16+ Non-Classical Monocytes"
                     "PPBP", "PF4", "ITGA2B", "GP9", "TUBB1", # "Platelets"
                     "FCER1A", "CD1C", "CLEC10A", "HLA.DRA", "ITGAX" # CD11c = ITGAX; "FCER1A+ Dendritic cells"
)

# PI3K/Akt pathway  
# c("IRS1", "PIK3CD", "PIK3R1", "PIK3AP1", "PIK3C3", "PTEN", 
#                "MTOR", "RICTOR", "PDK1", "AKT1", "AKT2", "AKT3", "MDM2", "TP53", 
#                "TSC2", "RHEB", "GSK3A", "GSK3B", "FOXO3", "TBK1", "RPTOR", 
#                "TRAF6", "NFKB1", "RELA", "REL", "RELB")

# TCA cycle genes
# c("PDHA1", "ACO2", "IDH1", "OGDH", "SUCLA2", "SUCLG1", "SDHA", "SDHC", "FH", "MDH1", "CS")

# ETC genes
# c("CYB561", "ENOX2", "ETFA", "ETFB", "ETFDH", "ETFRF1", 
#                "NDOR1", "NDUFA5", "PDIA5", "POR", "SDHA", "SDHB", "UQCRFS1", 
#                "DHRS2", "PPARGC1A")

PATH_output_vln <- paste0(PATH_output_figures, "/specific_mrkr_vln/")
dir.create(PATH_output_vln, recursive = TRUE)
for(ftr in mrkrs.2.plot){
  print(ftr)
  
  p <- VlnPlot(seurat.object, features = ftr,
               pt.size = 0.1) +         
    # geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +  # Add box plot manually
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 90, size = 6))
  print(p)
  ggsave(filename = paste0(PATH_output_vln, ftr, "_vln_plot.pdf"),
         plot = p, width = 12, height = 8, device = 'pdf', units = "in", useDingbats = FALSE)
}

########
# Create bubble plots
########
dir.create(paste0(PATH_output_figures, "/bubble_plots"))

# Cluster bubble plot
# Count number of cells per expID and cluster
cell_counts <- seurat.object@meta.data %>%
  group_by(expID, cluster = Idents(seurat.object)) %>%
  summarise(count = n(), .groups = "drop")

# Calculate proportions per expID
cell_counts <- cell_counts %>%
  group_by(expID) %>%
  mutate(proportion = count / sum(count)) # Normalize within expID

# Step 1: Compute the dominant cluster proprotion per expID
expID_order <- cell_counts %>%
  group_by(expID) %>%
  slice_max(order_by = proportion, n = 1) %>% # Select the cluster with the highest proportion per expID
  arrange(desc(proportion)) %>% # Sort by highest proportion
  pull(expID)

# Ensure clusters are ordered
# cell_counts$cluster <- factor(cell_counts$cluster, levels = unique(cell_counts$cluster))
# cell_counts$expID <- factor(cell_counts$expID, levels = expID_order)

# Bubble plot
p <- ggplot(cell_counts, aes(x = cluster, y = expID, size = proportion, fill = factor(cluster))) +
  geom_point(shape = 21, color = "black", alpha = 0.8) + # Bubble outline and transparency
  scale_size_area(max_size = 10) + # Adjust max bubble size
  scale_fill_viridis_d(option = "plasma") + # Color scale
  labs(x = "Cluster", y = "Experiment ID", size = "Proportion") +
  theme(
    legend.position = "none",
    # legend.box = "vertical",
    legend.margin = margin(10, 10, 10, 10),  # Add padding to prevent cropping
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate cluster labels
    panel.grind.major = element_line(color = "grey80"),
    panel.grind.minor = element_blank()
  )

print(p)
ggsave(filename = paste0(PATH_output_figures, "/bubble_plots/AML_bubble_plot.pdf"), 
       plot = p, width = 12, height = 8, device = 'pdf', units = "in", useDingbats = FALSE)

# Create top DEG bubble plot
ngenes <- 3 # Set number of top DEGs to look at

# Select top ngenes DEGs for each cluster
top.cluster.DEGs <- l.new.markers[[1]] %>% group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n=ngenes) %>% ungroup() %>%
  arrange(cluster, desc(avg_log2FC)) %>% # Sort genes within each cluster by logFC
  mutate(gene= factor(gene, levels = unique(gene))) # Convert to ordered factor

p <- ggplot(top.cluster.DEGs, aes(x = cluster, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21, color = "black", alpha = 0.8) +  # Bubble outline and transparency
  scale_size_area(max_size = 10) +  # Bubble size based on pct.1
  scale_fill_viridis_c(option = "plasma") + # Color scale
  labs(x = "Cluster", y = "Top DE Genes", size = "Proportion Expressed", fill = "avg_log2FC") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate cluster labels
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank()
  )

print(p)
ggsave(filename = paste0(PATH_output_figures, "/bubble_plots/gene_bubble_plot.pdf"), 
       plot = p, width = 12, height = 5*ngenes, device = 'pdf', units = "in", useDingbats = FALSE)



# Plot canonical marker bubble plot per cluster
# mrkrs.2.plot defined in the violin plot section
path <- paste0(PATH_output_figures, "/bubble_plots/canonical_mrkrs.pdf")

bubble.data <- l.new.markers[[1]] %>%
  filter(gene %in% mrkrs.2.plot) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  mutate(cluster = factor(cluster),
         gene = factor(gene, levels = mrkrs.2.plot))

p <- ggplot(bubble.data, aes(x = cluster, y = gene, size = pct.1, fill = avg_log2FC)) +
  geom_point(shape = 21, color = "black", alpha = 0.8) +  # Bubble outline and transparency
  scale_size_area(max_size = 10) +  # Bubble size based on pct.1
  scale_fill_viridis_c(option = "plasma") + # Color scale
  labs(x = "Cluster", y = "Markers", size = "Proportion Expressed", fill = "avg_log2FC") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate cluster labels
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank()
  )

print(p)
if(length(mrkrs.2.plot) < 6){
  ggsave(filename = path, plot = p,
         width = 12, height = 2*length(mrkrs.2.plot), 
         device = 'pdf', units = "in", useDingbats = FALSE)
} else {
  ggsave(filename = path, plot = p,
         width = 12, height = length(mrkrs.2.plot),
         device = 'pdf', units = "in", useDingbats = FALSE)
}

#########
# Volcano plot
#########
for(clst in unique(l.new.markers[[1]]$cluster)){
  mrkrs <- l.new.markers[[1]] %>% filter(cluster == clst)
  
  # The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
  mrkrs$diffexpressed <- "NO" # add a column of NAs
  
  mrkrs$diffexpressed[mrkrs$avg_log2FC > 0.5 & mrkrs$p_val_adj < 0.05] <- "UP" # If avg_log2FC > 0.6 and p_val < 0.05, set as "UP"
  mrkrs$diffexpressed[mrkrs$avg_log2FC < -0.5 & mrkrs$p_val_adj < 0.05] <- "DOWN" # If avg_log2FC < -0.6 and p_val < 0.05, set as "DOWN"
  
  mrkrs$label <- NA # Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
  mrkrs$label[mrkrs$diffexpressed != "NO"] <- mrkrs$gene[mrkrs$diffexpressed != "NO"]
  
  volcano.colors <- c("blue", "red", "black")
  names(volcano.colors) <- c("DOWN", "UP", "NO")
  
  # Subset up- and down-regulated genes
  up_genes <- mrkrs %>% 
    filter(diffexpressed == "UP") %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 10)
  
  down_genes <- mrkrs %>%
    filter(diffexpressed == "DOWN") %>%
    arrange(avg_log2FC) %>%
    slice_head(n = 10)
  
  # Combine top up- and down-regulated genes
  top_labels <- bind_rows(up_genes, down_genes)
  
  mrkrs$label <- NA
  mrkrs$label[mrkrs$gene %in% top_labels$gene] <- mrkrs$gene[mrkrs$gene %in% top_labels$gene]
  
  # Replace zero p-values with a tiny number before plotting
  mrkrs$p_val_adj[mrkrs$p_val_adj == 0] <- 1e-300
  
  x_max <- max(abs(mrkrs$avg_log2FC), na.rm = TRUE)
  y_max <- max(-log10(mrkrs$p_val_adj), na.rm = TRUE)
  
  # Add extra space so labels don't run off the edge
  x_pad <- 0.1 * x_max
  y_pad <- 0.1 * y_max
  
  p <- ggplot(data=mrkrs, aes(x=avg_log2FC, y=-log10(p_val_adj), 
                              col=diffexpressed, label=label)) + 
    geom_point() + 
    geom_vline(xintercept = c(-0.5,0.5), col="red") +
    geom_hline(yintercept = -log10(0.05), col="red") + 
    scale_color_manual(values=volcano.colors) + 
    geom_text_repel(
      na.rm = TRUE,
      color = "black",
      segment.color = "black",
      max.overlaps = Inf,
      box.padding = 0.75,        # Increase padding around labels
      point.padding = 0.5,       # Padding around points
      force = 2,                 # Repulsion strength
      force_pull = 0.5,          # Force pulling labels toward their point
      min.segment.length = 0.2   # Minimum length for connectors
    ) +
    coord_cartesian(xlim = c(-x_max - x_pad, x_max + x_pad),
                    ylim = c(0, y_max + y_pad))
    p <- p + expand_limits(x = c(-x_max - x_pad, x_max + x_pad), 
                           y = c(0, y_max + y_pad)) +
      ggtitle(paste(clst, "DEG volcano plot")) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")  # Center and bold the title
      )
  clst_name <- gsub(" ", "_", clst)
  ggsave(filename = paste0(PATH_output_figures, "/volcano_plots/", clst_name,
                           "_volcano.pdf"), plot = p, device = 'pdf', units = "in", useDingbats = FALSE)
}

#########
# Gene Set Enrichment Analysis
#########
l.degs <- split(l.new.markers[[1]], l.new.markers[[1]]$cluster)

for(cl in unique(l.new.markers[[1]]$cluster)){
  print(cl)
  deg.df <- subset(l.new.markers[[1]], cluster = cl)
  print(unique(deg.df$cluster))
}
ranks <- l.new.markers[[1]]$avg_log2FC
names(ranks) <- l.new.markers[[1]]$gene

#########
# Gene ontology
#########
GO_path = paste0(PATH_output_figures, "/GO_plots")

# Subselect clusters if needed
filt.mrkrs.4.go <- l.new.markers[[2]] %>%
  filter(grepl("NPM", cluster, ignore.case = FALSE))

l_GO_results <- calculate_GO_terms(l.new.markers[[1]], species, 3, GO_path)
saveRDS(l_GO_results, file = paste0(PATH_output_tables, "/l_GO_results.Rds"))

# l_GO_results <- readRDS(paste0(PATH_output_tables, "/l_GO_results.Rds"))

# Loop through each cluster in List_allRes and create plots
plots <- lapply(names(l_GO_results), function(cluster_name) {
  cluster_results <- l_GO_results[[cluster_name]]
  plot_go_terms_for_cluster(cluster_results, cluster_name, 20, GO_path)
})

# Optionally, print all plots (this will display them in the console if running interactively)
lapply(plots, print)

print(paste("GO plots saved to: ", GO_path))

"Section to be corrected"

# Select specific terms in the GO pathways to be looked (can be expanded as needed)
target_terms <- c("mitochondria", "ER", "endoplasmic reticulum", "nucleus", "nuclear",
                  "glycolysis", "glycolytic", "gluconeogenesis", "acetyl-CoA", "NAD",
                  "NADH", "NADP", "NADPH", "pyruvate")

# Function to filter GO terms related to target pathways
filter_target_terms <- function(go_results){
  go_results[grep(paste(target_terms, collapse = "|"), go_results$Term, ignore.case = TRUE), ]
}

GO_target_plots <- lapply(names(l_GO_results), function(cluster_name){
  cluster_results <- l_GO_results[[cluster_name]]
  
  # Filter for target-related terms
  target_results <- filter_target_terms(cluster_results)
  
  # If there are target-related results, plot them
  if(nrow(target_results) > 0) {
    plot_go_terms_for_cluster(target_results, cluster_name, 20, paste0(PATH_output_figures, "/GO_plots/glycolysis"))
  } else {
    message(paste("No target-related GO terms found for cluster:", cluster_name))
  }
})
lapply(GO_target_plots, print)

##########
# Create annotation dataframe
#########
meta_data <- seurat.object@meta.data #Select seurat object meta.data
meta_data$cell_type <- Idents(seurat.object) #Add cluster names as a column

for(exp_id in unique(meta_data$expID)){
  print(exp_id)
  path <- paste0(PATH_input_IRIS_imaging, "/", exp_id, "/metadata/seurat_metadata.csv")
  md <- meta_data[meta_data$expID == exp_id,]
  write.csv(md, file = path, row.names = FALSE)
}

write.csv(meta_data, file = paste0(PATH_output_tables, "/seurat_metadata.csv"), row.names = FALSE)
print(paste("Final annotation dataframe saved to: ", 
            paste0(PATH_output_tables, "/seurat_meta.csv")))

##########
# Create signature dataframes and lists
##########
l.new.markers <- readRDS(file = paste0(PATH_output_tables, "/l_markers.rds"))

#Order and filter cluster DEGs
l.new.markers[[2]] %>%
  group_by(cluster) %>%
  filter(avg_log2FC >= 0.2) %>%
  top_n(n = 500, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) -> new.markers.ordered

# Subselect clusters if needed
filtered.mrkrs <- new.markers.ordered %>%
  filter(grepl("health", cluster, ignore.case = TRUE))

write.table(new.markers.ordered, file = paste0(PATH_output_tables, "/healthy_sig_top500.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


##########################
# Code is useless from here on out
##########################
########
# Plot refined signatures as z-scores for verification/validation
########
"Use AddModuleScore"
seurat.object <- readRDS(file = paste0(PATH_output_objects, "/seurat_object_zscore.Rds"))

#Score each cell for cell type signatures
# No need to re-run this segment if multiple scores being calculated, skip to Calculate and z-score
RNA.norm <- as.data.frame(GetAssayData(seurat.object, assay = "RNA"))
scaled <- apply(RNA.norm, 1, scale)
rownames(scaled) <- colnames(RNA.norm)
scaled <- t(scaled)

# Load in target signature
target.sig <- read.delim("/home/tferrari/updepla/users/tferrari/Experiments/RNA_Analysis/IRIS/AML_experiments/signatures/npm1_signatures/refined_signatures/npm1_cohort_sig.txt", 
                         sep='\t', header = TRUE)

# Calculate and store z-score
scaled.sig <- subset(scaled, rownames(scaled) %in% target.sig$gene)

#Check for NAs and switch their values to 0 if necessary
scaled.sig.per.cell <- colSums(scaled.sig, na.rm = TRUE) / nrow(scaled.sig)
scaled.sig.per.cell[is.nan(scaled.sig.per.cell)] <- 0

seurat.object@meta.data[, "npm1_cohort_sig"] <- scaled.sig.per.cell

plots <- FeaturePlot(seurat.object, features = "npm1_cohort_sig",
                     reduction = "umap", cols = c("lightgrey", "deepskyblue", "blue"),
                     min.cutoff = c("q25"), pt.size = 0.5, order = TRUE) +
  theme(
    plot.title = element_blank(),
    plot.margin = margin(20, 20, 30, 20)
  )
p <- wrap_plots(plots)
print(p)
ggsave(filename = paste0(PATH_output_figures, "/z-score_plots/npm1_cohort_q25.pdf"), 
       plot = p, width = 8, height = 8, device = 'pdf', units = "in", useDingbats = FALSE)
########
# Create exhaustive signature marker list for future z-score annotations
########
# Load in specific signatures
sig_txt_files <- list.files("/home/tferrari/updepla/users/tferrari/Experiments/RNA_Analysis/IRIS/AML_experiments/signatures/individual_signatures", 
                            pattern = "\\.txt$", full.names = TRUE)
merged_df <- NULL # Initialize empty df for merging of signature df

for(file in sig_txt_files) {
  if(grepl("\\.txt$", file)) { # Ensure file has ".txt" extension
    message(paste("Loading: ", file))
    
    data <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
    
    # Merge/Concatenate the dataframes
    if(is.null(merged_df)){
      merged_df <- data # First file initializes the merged dataframe
    } else {
      merged_df <- bind_rows(merged_df, data)
    }
  }
}
print(dim(merged_df))
head(merged_df)

write.table(merged_df, file = paste0("/home/tferrari/updepla/users/tferrari/Experiments/RNA_Analysis/IRIS/AML_experiments/signatures/merged_sig_df.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# 1. Identify unique clusters for each dataframe

unique_clusters_jp299 <- jp299.markers %>% filter(!cluster %in% ng064.markers$cluster)
unique_clusters_ng064 <- ng064.markers %>% filter(!cluster %in% jp299.markers$cluster)

# 2. Identify overlapping clusters and intersect their genes
overlapping_clusters <- intersect(jp299.markers$cluster, ng064.markers$cluster)

# For overlapping clusters, find intersecting genes
intersecting_rows <- jp299.markers %>%
  filter(cluster %in% overlapping_clusters) %>%
  inner_join(ng064.markers, by = c("cluster", "gene"))

# Average the columns with the same base name
intersecting_rows_averaged <- intersecting_rows %>%
  mutate(
    p_val = (p_val.x + p_val.y) / 2,
    avg_log2FC = (avg_log2FC.x + avg_log2FC.y) / 2,
    pct.1 = (pct.1.x + pct.1.y) / 2,
    pct.2 = (pct.2.x + pct.2.y) / 2,
    p_val_adj = (p_val_adj.x + p_val_adj.y) / 2
  ) %>%
  dplyr::select(cluster, gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)  # Select only averaged columns

# 3. Combine the unique rows and intersecting rows into a new dataframe
result <- bind_rows(unique_clusters_jp299, unique_clusters_ng064, intersecting_rows_averaged)
write.table(result, file = "/home/tferrari/updepla/users/tferrari/Experiments/RNA_Analysis/IRIS/AML_experiments/extended_signature_df.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

create_signature_list <- function(df) {
  #Filter rows where avg_log2FC >= 0.5
  filtered_df <- df 
  # Split the filtered dataframe into a list of dataframes by cluster
  cluster_list <- split(filtered_df, filtered_df$cluster)
  
  return(cluster_list)
}
# Save the extended and intersected signatures into a list as an "*.Rds" file
extended.PBMC.sig.l <- create_signature_list(result)
saveRDS(extended.PBMC.sig.l, 
        file = "/home/tferrari/updepla/users/tferrari/Experiments/RNA_Analysis/IRIS/AML_experiments/extended_signature_list.Rds")

# Print specific colors on the UMAP
cluster_colors <- rep("black", 16) # Initialize all clusters to black
cluster_colors[7] <- "red"         # Cluster 7 in red
cluster_colors[c(3, 11)] <- "gray" # Clusters 3 and 11 in gray
cluster_colors[c(5, 6, 14)] <- "orange" # Clusters 5, 6, and 14 in orange
# Assuming your Seurat object is named 'seurat_object'
p <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", cols = cluster_colors) +
  scale_color_manual(values = cluster_colors) + NoLegend()
print(p)
ggsave(filename = paste0(PATH_output_figures, "/colored_umap.pdf"), 
       plot = p, width = 8, height = 8, device = 'pdf', units = "in", useDingbats = FALSE)
