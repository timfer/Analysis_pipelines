#Author: Joern Pezoldt
#Date init: 12.01.2024
#Scopes:
## 1) Embed all CD8 ER data processed with IRIS
## 2) Superimposte image data classification for DEG analysis and visualization

#libraries
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
sample_ids_IRIS <- c("NG064")
sample_ids <- c(sample_ids_IRIS)
#Define paths
PATH_input_IRIS_sequencing <- paste0("/home/",userID,"/updepla/projects/iris/4_sequencing_analysis")
PATH_input_IRIS_imaging <- paste0("/home/",userID,"/updepla/projects/iris/6_imaging_analysis")

PATH_output <- paste0("/home/",userID,
                      "/updepla/users/", userID, 
                      "/Experiments/RNA_analysis/IRIS/CHUV_samples/IRIS_open_source/NG064_NPM1")
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
# signatures.PBMCs <- read.delim(paste0(PATH_signature,
#                                       "/AML_v_HD_cell_sig_JP299.txt"),
#                                sep = "\t")
signatures.PBMCs <- read.delim("/home/tferrari/updepla/users/tferrari/Experiments/RNA_Analysis/IRIS/CHUV_samples/signatures/top100_combined_sig.txt",
                               sep = "\t", header = TRUE)
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
split.layers <- "yes"
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
#Load open source count matrices
#####
PATH_open_source_mtx <- "/home/tferrari/updepla/users/tferrari/Experiments/RNA_Analysis/open_source/GSE116256_van_Galen/count_matrix"
data <- as.data.frame(fread(paste0(PATH_open_source_mtx, "/GSM3587940_AML329-D0.dem.txt.gz")))

# Convert data frame to matrix
rownames(data) <- data[[1]] # Set first column as row names (gene names)
data <- data[, -1] # Remove the gene column as they are now rownames
count_matrix <- as.matrix(data) # Create count matrix for the target dataset

cell_barcodes <- colnames(count_matrix) # Extract cell barcodes (column names of count matrix)
# Create a data frame for metadata
meta.data <- data.frame(
  Full_Cell_ID = cell_barcodes,  # Store the cell barcodes for reference
  Sample_ID = "AML210A-D0",  # Assign sample ID (modify as needed)
  n_UMI = colSums(count_matrix),  # Total counts per cell
  n_genes = colSums(count_matrix > 0),  # Number of detected genes per cell
  row.names = cell_barcodes  # Set row names to cell barcodes
)

# Combine the open source & IRIS meta.data
all_columns <- union(colnames(l.meta.mtx[[1]]), colnames(meta.data)) # Step 1: Identify all unique column names
# Step 2: Add missing columns to the open source meta.data, filling with NAs
for(col in setdiff(all_columns, colnames(meta.data))){
  meta.data[[col]] <- NA
}

# Step 3: Ensure the same column order in both dataframes
meta.data <- meta.data[, all_columns]
l.meta.mtx[[1]] <- l.meta.mtx[[1]][, all_columns]

combined_meta_data <- rbind(meta.data, l.meta.mtx[[1]])

# Combine the open source & IRIS count matrices
rownames(count_matrix) <- gsub("^MT\\.", "MT-", rownames(count_matrix))
rownames(l.meta.mtx[[2]]) <- gsub("^MT\\.", "MT-", rownames(l.meta.mtx[[2]]))

all_genes <- union(rownames(count_matrix), rownames(l.meta.mtx[[2]])) # Step 2: Identify all unique genes
expand_mtx <- function(mat, all_genes) { # Step 3: Expand matrices to include all genes
  missing_genes <- setdiff(all_genes, rownames(mat))
  if(length(missing_genes) > 0) {
    # Create a matrix of zeros for missing genes
    zero_mtx <- matrix(0, nrow = length(missing_genes), ncol = ncol(mat))
    rownames(zero_mtx) <- missing_genes
    colnames(zero_mtx) <- colnames(mat)
    # Combine the original matrix with the zero matrix
    mat <- rbind(mat, zero_mtx)
  }
  # Reorder rows to match all_genes
  mat <- mat[all_genes, ]
  return(mat)
}
# Expand both matrices
count_mtx_expanded <- expand_mtx(count_matrix, all_genes)
mtx.expanded <- expand_mtx(l.meta.mtx[[2]], all_genes)

# Step 4: Combine the matrices
combined.mtx <- cbind(count_mtx_expanded, mtx.expanded)

#####
#Batch integration with Seurat
#####
#Seurat
integration.results <- seurat_integrate(
  signature_cell_cycle_human_mouse,
  species,
  cond.2.rmv,
  combined.mtx,
  combined_meta_data,
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
dims_use <- integration.results[[2]]
integration.name <- integration.results[[3]]

resolutions <- seq(0.5, 1.0, by = 0.1) # Define the range of resolutions
find_cluster_resolution(seurat.object, dims_use, resolutions, integration.name,
                        cluster_algo, paste0(PATH_output_figures, "/UMAP_plots"))

if(!file.exists(paste0(PATH_output_figures, "/unnamed_plots"))){
  dir.create(paste0(PATH_output_figures, "/unnamed_plots"), recursive = TRUE)
  message(paste("Directory created: ", paste0(PATH_output_figures, "/unnamed_plots")))
} else {
  message(paste("Directory already exists: ", paste0(PATH_output_figures, "/unnamed_plots")))
}

seurat.object <- cluster_UMAP_seurat(seurat.object, 0.7, dims_use, integration.name,
                                     reduction.name, cluster_algo, 
                                     paste0(PATH_output_figures, "/unnamed_plots"))
#Save seurat.object
saveRDS(file = paste0(PATH_output_objects, "/seurat_object_integrated.Rds"),
        seurat.object)
print(paste("Fully integrated seurat object saved to: ",
            paste0(PATH_output_objects, "/seurat_object_integrated.Rds")))

# seurat.object <- readRDS(file = paste0(PATH_output_objects, "/seurat_object_integrated.Rds"))

#############################
#Plotting
#############################
#####
#QC
#####
generate_QC_plots(seurat.object, reduction.name, c("NG064"), 
                  paste0(PATH_output_figures, "/QC_plots"))
print(paste("QC plots saved to: ", paste0(PATH_output_figures, "/QC_plots")))

generate_bar_plots(seurat.object, paste0(PATH_output_figures, "/QC_plots"))
print(paste("Barplots saved to: ", paste0(PATH_output_figures, "/QC_plots")))

#####
#UMAP
#####
generate_ftr_plots_umap(seurat.object, reduction.name,
                        paste0(PATH_output_figures, "/unnamed_plots"))
print(paste("Feature and UMAP plots saved to: ",
            paste0(PATH_output_figures, "/unnamed_plots")))

#####
#PBMC signatures
#####
#Select top 500 DEGs per cluster and order in descending for avg_log2FC and plot
seurat.object <- JoinLayers(seurat.object)

signatures.PBMCs %>% group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.25) %>%
  slice_head(n = 250) %>%
  ungroup() -> top.PBMC.DEGs

gene_list <- list()
for(score_name in unique(signatures.PBMCs$cluster)){
  print(score_name)
  gene_list[[score_name]] <- top.PBMC.DEGs[top.PBMC.DEGs$cluster == score_name, ]$gene
  seurat.object <- AddModuleScore(seurat.object, 
                                  features = list(gene_list[[score_name]]), 
                                  name = score_name)
}

# Remove the "1" suffix from the new score columns in meta.data
# Renaming columns in the meta.data slot
names(seurat.object@meta.data)[(ncol(seurat.object@meta.data) - length(gene_list) + 1):ncol(seurat.object@meta.data)] <- names(gene_list)

p <- FeaturePlot(seurat.object, features = names(gene_list),
                 reduction = "umap", cols = c("lightgrey", "deepskyblue", "blue"),
                 min.cutoff = "q70", pt.size = 0.5, order = TRUE)
print(p)
ggsave(filename = paste0(PATH_output_figures, 
                         "/unnamed_plots",
                         "/cell_type_ftr_plot.png"), plot = p, width = 24, height = 32)

seurat.object <- calculate_cell_signatures(signatures.PBMCs, seurat.object, 100,
                                           reduction.name, paste0(PATH_output_figures, 
                                                                  "/unnamed_plots",
                                                                  "/cell_type_ftr_plot.png"))

saveRDS(file = paste0(PATH_output_objects, "/seurat_object_zscore.Rds"), seurat.object)
print(paste("seurat object with Z-scores in @meta_data saved to: ",
            paste0(PATH_output_objects, "/seurat_object_zscore.Rds")))

# seurat.object <- readRDS(file = paste0(PATH_output_objects, "/seurat_object_zscore.Rds"))

################################################
#DEGs
################################################
#####
#DEGs general
#####
l.markers <- calculate_and_plot_marker_DEGs(seurat.object, min.pct,
                                            logfc.threshold, 1000,
                                            paste0(PATH_output_figures,
                                                   "/unnamed_plots/DEG_heatmap.png"))
# Also save it to an excel file for easy reading
wb <- createWorkbook() #Create new workbook
addWorksheet(wb, "Unnamed_all_markers") #Add worksheet
writeData(wb, "Unnamed_all_markers", l.markers[[1]]) #Write data
addWorksheet(wb, "Unnamed_top_markers") #Add worksheet
writeData(wb, "Unnamed_top_markers", l.markers[[2]]) #Write data
addWorksheet(wb, "Unnamed_bottom_markers") #Add worksheet
writeData(wb, "Unnamed_bottom_markers", l.markers[[3]]) #Write data
saveWorkbook(wb, paste0(PATH_output_tables, "/unnamed_markers.xlsx")) #Save workbook

########
# Feature plots
########
plot_specific_ftr(seurat.object, paste0(PATH_output_figures, "/feature_plots"))
print(paste("Class specific feature plots saved to ", 
            paste0(PATH_output_figures, "/feature_plots")))

########
# Violin plots for specific blast/oncological markers
########
#Subset the seurat object to contain only the targeted experiment
PATH_output_vln <- paste0(PATH_output_figures, "/vln_plots")
dir.create(PATH_output_vln, recursive = TRUE)

#Find genes with specific sequences
# Extract all gene names from the Seurat object
all_genes <- rownames(seurat.object)
target_genes <- grep("ATP", all_genes, value = TRUE, ignore.case = TRUE)

ftrs.2.plot <- c("WT1")

# TCA cycle genes
# c("PDHA1", "ACO2", "IDH1", "OGDH", "SUCLA2", "SUCLG1", "SDHA", "SDHC", "FH", "MDH1", "CS")

# ETC genes
# c("CYB561", "ENOX2", "ETFA", "ETFB", "ETFDH", "ETFRF1", 
#                "NDOR1", "NDUFA5", "PDIA5", "POR", "SDHA", "SDHB", "UQCRFS1", 
#                "DHRS2", "PPARGC1A")

for(ftr in ftrs.2.plot){
  p <- VlnPlot(seurat.object, features = ftr) + theme(legend.position = 'none')
  print(p)
  ggsave(filename = paste0(PATH_output_vln, ftr, "_vln_plot.jpeg"),
         plot = p, width = 12, height = 8)
}

######
#Cluster annotation
######
#Rename clusters with major cell group attribution and then subselect again
#Do not use the same names as for signature_PBMCs

# NPM1 v HD cluster annotation
new.cluster.names <-  c("CD14 Classical Monocytes", "NPM1 AML blasts", # Cluster 1-2
                        "CD56 NK Cells", "CD4 Central Memory T Cells", # Cluster 3-4
                        "Effector Th2 Cells", "CD8 Cytotoxic T Cells", # Cluster 5-6
                        "NPM1 Differentiating Blasts", "Naive CD4 T Cells", # Cluster 7-8
                        "CD8 Effector Memory T Cells", "Naive CD8 T Cells", # Cluster 9-10
                        "Healthy CD19 B Cells", "CD16 Non-Classical Monocytes", # Cluster 11-12
                        "Dendritic Cells", "Activated NK Cells") # Cluster 13-14
# AML-MRC v HD cluster annotation 
# c("CD14 Classical Monocytes", "CD56 NK Cells", # Cluster 1-2
#                      "CD4 Central Memory T Cells", "CD8 Cytotoxic T Cells", # Cluster 3-4
#                      "Effector Th2 Cells", "AML-MRC LSPCs", # Cluster 5-6
#                      "Naive CD8 T Cells", "Naive CD4 T Cells", # Cluster 7-8
#                      "AML-MRC Erythroblast", "CD8 Effector Memory T Cells", # Cluster 9-10
#                      "AML-MRC CD19 B Cells", "Healthy CD19 B Cells", # Cluster 11-12
#                      "CD16 Non-Classical Monocytes", "Dendritic Cells", # Cluster 13-14
#                      "AML-MRC Proliferating/Cycling Blasts", "Activated NK Cells", # Cluster 15-16
#                      "AML-MRC-associated Monocytes") # Cluster 17

# NG064 v JP299 v HD
# c("CD14 Classical Monocytes", "NPM1 LSPCs", # Cluster 1-2
#   "CD56 NK cells","CD4 Central Memory T Cells", # Cluster 3-4
#   "Effector Th2 Cells","CD8 Cytotoxic T Cells", # Cluster 5-6
#   "AML-MRC LSPCs", "Naive CD8 T Cells", # Cluster 7-8
#   "NPM1 Differentiating Myeloid Progenitors","CD8 Effector Memory T Cells", # Cluster 9-10
#   "AML-MRC Erythroblasts","Naive CD4 T Cells with Early Activation", # Cluster 11-12
#   "AML-MRC CD19 B Cells","Healthy CD19 B Cells", # Cluster 13-14
#   "CD16 Non-Classical Monocytes","Dendritic Cells", # Cluster 15-16
#   "AML-MRC Proliferating Blasts","Activated NK Cells", # Cluster 17-18
#   "AML-MRC-associated Monocytes") # Cluster 19

seurat.object <- rename_clusters(seurat.object, new.cluster.names, reduction.name,
                                 paste0(PATH_output_figures, "/named_plots"))

saveRDS(file = paste0(PATH_output_objects, "/seurat_object_final.Rds"), seurat.object)
print(paste("Final fully integrated seurat object saved to: ",
            paste0(PATH_output_objects, "/seurat_object_final.Rds")))
# seurat.object <- readRDS(file = paste0(PATH_output_objects, "/seurat_object_final.Rds"))

##########
# Detect new markers
##########
l.new.markers <- calculate_and_plot_marker_DEGs(seurat.object, min.pct, logfc.threshold, 
                                                1000, paste0(PATH_output_figures,
                                                             "/named_plots/named_DEG_heatmap.png"))
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

saveRDS(l.new.markers, file = paste0(PATH_output_tables, "/l_markers.Rds"))
print(paste("Markers list saved to: ", paste0(PATH_output_tables, "/l_markers.Rds")))
print("Note: Markers list 1: All markers, 2: top 500 EDGs, 3: bottom 500 DEGs")

# l.new.markers <- readRDS(file = paste0(PATH_output_tables, "/l_markers.Rds"))

#########
# Volcano plot
#########
# Subset for target clusters
aml.mrkrs <- l.new.markers[[1]] %>% filter(cluster %in% c("NPM1 LSPCs", "NPM1 Differentiating Myeloid Progenitors"))
volcano_plot <- ggplot(aml.mrkrs, aes(x = avg_log2FC, y = -log(p_val),
                                      color = ifelse(p_val_adj > -log(0.1) & abs(avg_log2FC) > 2,
                                                     "red", "black"))) +
  geom_point(size = 2) +
  labs(title = "Volcano Plot: NPM1-mut LSPCs v Differentiating Myeloid Cells",
       x = "Expression (avg_log2FC)", y = "-log10(p-value") +
  theme_minimal()

print(volcano_plot)
#########
# Gene ontology
#########
GO_path = paste0(PATH_output_figures, "/GO_plots")
l_GO_results <- calculate_GO_terms(l.new.markers[[2]], species, 3, GO_path)
saveRDS(l_GO_results, file = paste0(PATH_output_tables, "/l_GO_results.Rds"))

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
target_terms <- c("glycolysis", "glucose", "hexose", "fructose", "carbohydrate metabolism")

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
create_annotation_dataframe(seurat.object, l.image.doublet.ROI.annotation, 
                            paste0(PATH_output_tables, "/annotation_df.csv"))
print(paste("Final annotation dataframe saved to: ", 
            paste0(PATH_output_tables, "/annotation_df.csv")))

##########
# Create signature dataframes and lists
##########
#Order and filter cluster DEGs
l.new.markers[[2]] %>%
  group_by(cluster) %>%
  filter(avg_log2FC >= 0.2) %>%
  top_n(n = 500, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) -> new.markers.ordered

write.table(new.markers.ordered, file = paste0(PATH_output_tables, "/AML_v_HD_cell_sig_JP299.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

########
# Create exhaustive signature marker list for future z-score annotations
########
# 1. Identify unique clusters for each dataframe
jp299.markers <- read.table("/home/tferrari/updepla/users/tferrari/Experiments/RNA_Analysis/IRIS/CHUV_samples/JP299_a_CH2_AML/tables/AML_v_HD_cell_sig_JP299.txt",
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ng064.markers <- read.table("/home/tferrari/updepla/users/tferrari/Experiments/RNA_Analysis/IRIS/CHUV_samples/NG064_a_CH1_AML/tables/AML_v_HD_cell_sig_NG064.txt",
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)

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
write.table(result, file = "/home/tferrari/updepla/users/tferrari/Experiments/RNA_Analysis/IRIS/CHUV_samples/extended_signature_df.txt",
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
        file = "/home/tferrari/updepla/users/tferrari/Experiments/RNA_Analysis/IRIS/CHUV_samples/extended_signature_list.Rds")

########
# Select non-intersecting genes
########
# top.markers <- l.new.markers[[2]]
# aml.markers.list <- list(top.markers[top.markers$cluster == "Leukemic Stem and Progenitor Cells",], 
#                         top.markers[top.markers$cluster == "AML-MRC Erythroblast",], 
#                         top.markers[top.markers$cluster == "AML-MRC Proliferating/Cycling Blasts",])
# 
# aml.assoc.markers.l <- list(top.markers[top.markers$cluster == "AML-MRC CD19+ B Cells",],
#                             top.markers[top.markers$cluster == "AML-MRC-associated Monocytes",])
# 
# healthy.markers.l <- list(top.markers[top.markers$cluster == "CD14+ Classical Monocytes",],
#                           top.markers[top.markers$cluster == "CD56+ Natural Killer Cells",],
#                           top.markers[top.markers$cluster == "CD4+ Central Memory T Cells",],
#                           top.markers[top.markers$cluster == "CD8+ Cytotoxic T Cells",],
#                           top.markers[top.markers$cluster == "Effector Th2 Cells",],
#                           top.markers[top.markers$cluster == "Naive CD8+ T Cells",],
#                           top.markers[top.markers$cluster == "Naive CD4+ T Cells with Early Activation",],
#                           top.markers[top.markers$cluster == "CD8+ Effector Memory T Cells",],
#                           top.markers[top.markers$cluster == "Healthy CD19+ B Cells",],
#                           top.markers[top.markers$cluster == "CD16+ Non-Classical Monocytes",],
#                           top.markers[top.markers$cluster == "Dendritic Cells",],
#                           top.markers[top.markers$cluster == "Activated Natural Killer Cells",])
# 
# #Find all unique genes across all tibbles
# unique.genes <- aml.markers.list %>%
#   map(~ pull(., gene)) %>% # Extract the genes column from each tibble
#   reduce(union) # Find the union of all gene sets
# 
# # Find the intersected genes (those that appear in all tibbles)
# intersected.genes <- aml.markers.list %>%
#   map(~ pull(., gene)) %>% # Extract the genes column from each tibble
#   reduce(intersect) # Find the intersection of all gene sets
# 
# non.intersecting.genes <- setdiff(unique.genes, intersected.genes)
