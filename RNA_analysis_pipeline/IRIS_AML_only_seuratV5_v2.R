

library(dplyr)
library(patchwork)
library(Seurat)
library(SeuratData)
library(UpSetR)
library(grid)
library(gridExtra)
library(ggplot2)
library(topGO)
library(reticulate)
library(leidenAlg)
library(openxlsx)
library(RColorBrewer)
library(scales)

#load custom functions around scRNA-seq with IRIS, imaging and tools used for analysis
PATH_functions <- "/home/tferrari/NAS2/iris/1_scripts/iris_scRNAseq/2_functions/scRNA-seq/PBMC_seurat_functions"
source(paste0(PATH_functions,"/seuratV5_integration_v2.R"))
source("/home/tferrari/NAS2/iris/1_scripts/iris_scRNAseq/2_functions/scRNA-seq/support_transcriptome_integration.R")
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
sample_ids_IRIS <- c("JP299")
sample_ids <- c(sample_ids_IRIS)
#Define paths
PATH_input_IRIS_sequencing <- paste0("/home/",userID,"/updepla/projects/iris/4_sequencing_analysis")
PATH_input_IRIS_imaging <- paste0("/home/",userID,"/updepla/projects/iris/6_imaging_analysis")

PATH_output <- paste0("/home/",userID,
                      "/updepla/users/", userID, 
                      "/Experiments/RNA_analysis/IRIS/CHUV_samples/JP299/JP299_only")
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
PATH_signature <- "/home/tferrari/NAS2/iris/1_scripts/iris_scRNAseq/3_signatures_genes"
#Signatures PBMCs
signatures.PBMCs <- read.delim(paste0(PATH_signature,
                                      "/AML_PBMC_extended_markers.txt"),
                               sep = "\t")
# Signatures cell cycle
signature_cell_cycle_human_mouse <- readRDS(paste0(PATH_signature,
                                                   "/signature_cell_cycle_human_mouse.Rds"))

#####
#Parameters to be modified
#####
#Variables
species = "human"
# Batch over "CC" or "expID"
integrate_over = "expID"

#Threshing
min_cells_percent = 0.005
min_gene_number = 500
mito_cutoff = 15
min_nUMI = 3000

#Set starting dimensions
dims_use = 50
#common variables to regress are: "n_UMI", "percent.mt", "Phase"
vars_to_regress <- c()

set.seed(300)

#Set threshold values for marker DEGs
min.pct <- 0.1
logfc.threshold <- 0.2

#Set experiment IDs to be removed if necessary (only necessary if 
#want to maintain "uniqueness" of clusters for datasets with low relative cell counts)
cond.2.rmv <- c()

#Set whether to first split layers and then perform IntegrateLayers on the seurat.object
#Will split and integrate layers if "yes"
split.layers <- "no"
integrate.layers <- "no"
# Integration method for integrate layers: CCAIntegration, RPCAIntegration, HarmonyIntegration
# FastMNNIntegration, scVIIntegration
#integration.method <- CCAIntegration

#Select the type of community clustering algorithm
cluster_algo <- "Leiden"

#UMAP resolution
resolution <- 0.6
#Set UMAP reduction name
integration.name <- "pca"
reduction.name <- "umap"

# Prepare the text to be printed
text_content <- paste(
  "##### Parameters to be Modified #####\n",
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
  "Integration Method: ", deparse(substitute(integration.method)), "\n",
  "Clustering Algorithm: ", cluster_algo, "\n",
  "UMAP Resolution: ", resolution, "\n",
  "UMAP Integration Name: ", integration.name, "\n",
  "UMAP Reduction Name: ", reduction.name, "\n"
)

writeLines(text_content, con = paste0(PATH_output, "/integration_parameters.txt"))

#######################################
#Image data
#######################################
#####
#Load Count matrix IRIS platform
#####
#Annotate doublets and store them in a list
l.image.doublet.ROI.annotation <- annotate_doublets(PATH_input_IRIS_imaging, 
                                                    sample_ids_IRIS)

#LiquideDrop data
l.meta.mtx <- create_matrix_meta_data(sample_ids_IRIS, l.image.doublet.ROI.annotation)

#####
#Batch integration with Seurat
#####
#Seurat
print(unique(l.meta.mtx[[1]]$CellType_Condition))

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
  paste0(PATH_output_figures, "/QC_plots"))

seurat.object <- integration.results[[1]]

resolutions <- seq(1.0, 2.0, by = 0.1)# Define the range of resolutions
find_cluster_resolution(seurat.object, 38, resolutions,
                        paste0(PATH_output_figures, "/UMAP_plots"))

seurat.object <- cluster_UMAP_seurat(seurat.object, 1.1, 38, 
                                     reduction.name)
DimPlot(seurat.object)
#Save seurat.object
saveRDS(file = paste0(PATH_output_objects, "/seurat_object_integrated.Rds"), seurat.object)
print(paste("Fully integrated seurat object saved to: ",
            paste0(PATH_output_objects, "/seurat_object_integrated.Rds")))
seurat.object <- readRDS(paste0(PATH_output_objects, "/seurat_object_integrated.Rds"))
#############################
#Plotting
#############################
#####
#QC
#####
generate_QC_plots(seurat.object, reduction.name, "", paste0(PATH_output_figures, "/QC_plots"))
print(paste("QC plots saved to: ", paste0(PATH_output_figures, "/QC_plots")))
#####
#UMAP
#####
generate_ftr_plots_umap(seurat.object, reduction.name, paste0(PATH_output_figures,
                                                              "/unnamed_plots"))
print(paste("Feature and UMAP plots saved to: ", paste0(PATH_output_figures, "/unnamed_plots")))
#####
#PBMC signatures
#####
#Select top 500 DEGs per cluster and order in descending for avg_log2FC and plot
seurat.object <- calculate_cell_signatures(signatures.PBMCs, seurat.object, 
                                           paste0(PATH_output_figures, "/unnamed_plots",
                                                  "/cell_type_ftr_plot.png"))
seurat.object <- calculate_cell_signatures(signatures.AML, seurat.object,
                                           paste0(PATH_output_figures, "/unnamed_plots",
                                                  "/aml_ftr_plot.png"))

saveRDS(file = paste0(PATH_output_objects, "/seurat_object_zscore.Rds"), seurat.object)
print(paste("seurat object with Z-scores in @meta_data saved to: ",
            paste0(PATH_output_objects, "/seurat_object_zscore.Rds")))

################################################
#DEGs
################################################
#####
#DEGs general
#####
l.markers <- calculate_and_plot_marker_DEGs(seurat.object, min.pct, logfc.threshold, 
                                            paste0(PATH_output_figures,
                                                   "/unnamed_plots/DEG_heatmap.png"))
top500 <- l.markers[[2]]

#Also save it to an excel file for easy reading
wb <- createWorkbook() #Create new workbook
addWorksheet(wb, "Unnamed_top250_markers") #Add worksheet
writeData(wb, "Unnamed_top250_markers", top500) #Write data
saveWorkbook(wb, paste0(PATH_output_tables, "/unnamed_top250_markers.xlsx")) #Save workbook

########
# Feature plots
#######
plot_specific_ftr(seurat.object, paste0(PATH_output_figures, "/feature_plots"))
print(paste("Class specific feature plots saved to ", paste0(PATH_output_figures, "/feature_plots")))

######
#Cluster annotation
######
#Rename clusters with major cell group attribution and then subselect again
#Do not use the same names as for signature_PBMCs
new.cluster.names <- c("AML_stem_cells_1", "AML_blasts_1", "AML_healthy_memory_CD8_1", "AML_healthy_monocytes_1", # Cluster 0-3
                       "AML_erythroid_progenitors_1", "AML_healthy_effector_CD8_1", "AML_healthy_CD4_1", "AML_granulocytes_1", #Cluster 4-7
                       "AML_intermediate_diff_myeloid_1", "AML_healthy_B_cells_1") #Cluster 8-9
seurat.object <- rename_clusters(seurat.object, new.cluster.names, 
                                 paste0(PATH_output_figures, "/named_plots"))

saveRDS(file = paste0(PATH_output_objects, "/seurat_object_final.Rds"), seurat.object)
print(paste("Final fully integrated seurat object saved to: ",
            paste0(PATH_output_objects, "/seurat_object_final.Rds")))
#seurat.object <- readRDS(file = paste0(PATH_output_objects, "/seurat_object_final.Rds"))

##########
# Detect new markers
##########
l.new.markers <- calculate_and_plot_marker_DEGs(seurat.object, min.pct, logfc.threshold, 
                                                paste0(PATH_output_figures,
                                                       "/named_plots/named_DEG_heatmap.png"))
new.markers.all <- l.new.markers[[1]]
top500 <- l.new.markers[[2]]
bottom500 <- l.new.markers[[3]]

saveRDS(l.new.markers, file = paste0(PATH_output_tables, "/l_markers.Rds"))
print(paste("Markers list saved to: ", paste0(PATH_output_tables, "/l_markers.Rds")))
print("Note: Markers list 1: All markers, 2: top 500 EDGs, 3: bottom 500 DEGs")
##########
# Gene ontology
#########
GO_path = paste0(PATH_output_figures, "/GO_plots")
l_GO_results <- calculate_GO_terms(top500, species, GO_path)

saveRDS(l_GO_results, file = paste0(PATH_output_tables, "/l_GO_results.Rds"))

# Loop through each cluster in List_allRes and create plots
plots <- lapply(names(l_GO_results), function(cluster_name) {
  cluster_results <- l_GO_results[[cluster_name]]
  plot_go_terms_for_cluster(cluster_results, cluster_name, 20, GO_path)
})

# Optionally, print all plots (this will display them in the console if running interactively)
lapply(plots, print)

print(paste("GO plots saved to: ", GO_path))

##########
# Create annotation dataframe
#########
create_annotation_dataframe(seurat.object, l.image.doublet.ROI.annotation, 
                            paste0(PATH_output_tables, "/annotation_df.csv"))
print(paste("Final annotation dataframe saved to: ", 
            paste0(PATH_output_tables, "/annotation_df.csv")))

##########
# Create signature dataframes and lists
##########
#Order and filter cluster DEGs
new.markers.all %>%
  group_by(cluster) %>%
  filter(avg_log2FC >= 0.25) %>%
  top_n(n = 500, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) -> new.markers.ordered

remove_suffix <- function(df) {
  df$cluster <- gsub("_1$", "", df$cluster)
  return(df)
}
new.markers.ordered <- remove_suffix(new.markers.ordered)

#Save extended PBMC signature dataframe as .txt file
write.table(new.markers.ordered, file = paste0("/home/tferrari/updepla_storage/users/tferrari/Experiments/RNA_Analysis/PBMC_seuratv5/AML_only_20240731/results/tables", 
                                               "/AML_only_sig.txt"), 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

create_signature_list <- function(df) {
  #Filter rows where avg_log2FC >= 0.5
  filtered_df <- df 
  # Split the filtered dataframe into a list of dataframes by cluster
  cluster_list <- split(filtered_df, filtered_df$cluster)
  
  return(cluster_list)
}

extended.PBMC.sig.l <- create_signature_list(new.markers.ordered)