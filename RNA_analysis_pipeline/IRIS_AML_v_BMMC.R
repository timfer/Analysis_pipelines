#Author: Joern Pezoldt & T. Ferrari
#Date init: 12.01.2024
#Scopes:
## 1) Embed all CD8 ER data processed with IRIS
## 2) Superimposte image data classification for DEG analysis and visualization

"
Combined integrations even using seuratv5 split layer integrations, lead to forcing of the 
AML myeloblast populations into healthy BM clusters.
"

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
library(ggrepel)
library(EnhancedVolcano)
library(biomaRt)
library(data.table)

#load custom functions around scRNA-seq with IRIS, imaging and tools used for analysis
PATH_functions <- "/home/tferrari/updepla/users/tferrari/Experiments/RNA_Analysis/IRIS/scripts/PBMC_seurat_functions"
source(paste0(PATH_functions,"/seuratV5_integration_v2.R"))
source(paste0(PATH_functions, "/support_transcriptome_integration.R"))
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
sample_ids_IRIS <- c("NG078", "NG075",
                     "TF051"
                     )
sample_ids <- c(sample_ids_IRIS)
#Define paths
PATH_input_IRIS_sequencing <- paste0("/home/",userID,"/updepla/projects/iris/4_sequencing_analysis")
PATH_input_IRIS_imaging <- paste0("/home/",userID,"/updepla/users/tferrari/Experiments/Image_Analysis/AML_experiments/experiment_files")
PATH_gtf_annotation <- paste0("/home/", userID, "/updepla/users/", userID,
                              "/Experiments/RNA_Analysis/IRIS/gtf_annotations/gencode.v44.annotation.gtf")

PATH_output <- paste0("/home/",userID,
                      "/updepla/users/", userID, 
                      "/Experiments/RNA_analysis/IRIS/AML_experiments/results/BM-PB/everything/noosc_nolay")
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
signatures.PBMCs <- read.delim("/home/tferrari/updepla/users/tferrari/Experiments/RNA_Analysis/IRIS/AML_experiments/signatures/healthy_sig.txt",
                               sep = "\t", header = TRUE)
signatures.BMMCs <- read.delim("/home/tferrari/updepla/users/tferrari/Experiments/RNA_Analysis/IRIS/AML_experiments/signatures/NPM1_NG064_BMMC_TSP21_TSP25_sig.txt",
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
min_gene_number = 250
mito_cutoff = 15
min_nUMI = 2500

#Set starting dimensions
dims_use = 50
#common variables to regress are: "n_UMI", "percent.mt", "Phase"
vars_to_regress <- c()

set.seed(300)

#Set threshold values for marker DEGs
min.pct <- 0.05
logfc.threshold <- 0.25

#Set experiment IDs to be removed if necessary (only necessary if 
#want to maintain "uniqueness" of clusters for datasets with low relative cell counts)
cond.2.rmv <- c()

#Set whether to first split layers and then perform IntegrateLayers on the seurat.object
#Will split and integrate layers if "yes"
split.layers <- "yes"
integrate.layers <- "no"
integration.type <- "harmony"

#UMAP resolution
# resolution <- 0.8
#Set UMAP reduction name
integration.name <- "pca"
reduction.name <- "umap"

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

########
# Load in BMMC data and append the count matrix and meta_data to IRIS matrix and meta_data
########
seurat.BM <- readRDS(paste0("/home/", userID,
                            "/updepla/users/", userID,
                            "/Experiments/RNA_Analysis/IRIS/AML_experiments/open_source_controls/tabula_sapiens/",
                            "tab_sapiens_marrow.rds"))

# Create a count matrix and meta_data dataframe for combining with the IRIS data
mtx.BM <- as.matrix(seurat.BM@assays$RNA@counts)

# The seurat.BM rownames are in ENSG code and we need to convert it to HGNC
# Downloaded and unzipped gencode.v44.annotation.gtf.gz
gtf <- fread(PATH_gtf_annotation, skip = 5, sep = "\t", header = FALSE)
genes <- gtf[V3 == "gene", ] # Subset gtf to exclude exons and transcripts

# Extract gene_id and gene_name directly from V9
genes[, ensembl_id := sub('.*gene_id "([^"]+)".*', '\\1', V9)]
genes[, gene_type := sub('.*gene_type "([^"]+)".*', '\\1', V9)]
genes[, hgnc_id := sub('.*gene_name "([^"]+)".*', '\\1', V9)]

# Restrict according to gene type, e.g. exclude all pseudogenes etc.
print(unique(genes$gene_type))
patterns.2.exclude <- c("pseudogene", "lncRNA", "TEC", "misc_RNA", "artifact")
genes <- genes[!Reduce(`|`, lapply(patterns.2.exclude, function(p) grepl(p, genes$gene_type)))]

# 6. Remove version numbers from `ensembl_id`
genes[, ensembl_id := sub("\\..*", "", ensembl_id) ]

# 7. Drop duplicates
genes <- unique(genes[, .(ensembl_id, hgnc_id)])

# Create mapping
ensg2hgnc <- setNames(genes$hgnc_id, genes$ensembl_id)

# Match and rename rownames
new_gene_names <- ensg2hgnc[rownames(mtx.BM)]

# Filter out genes without HGNC symbols
valid <- !is.na(new_gene_names) & new_gene_names != ""
mtx.BM <- mtx.BM[valid, ]
rownames(mtx.BM) <- new_gene_names[valid]

# Optional: remove duplicated HGNC symbols (if present)
mtx.BM <- mtx.BM[!duplicated(rownames(mtx.BM)), ]

mtx <- l.meta.mtx[[2]] # IRIS count matrix
mtx <- mtx[!duplicated(rownames(mtx)), ]

# Create a count matrix as the union of all genes in each respective count matrix
all_genes <- union(rownames(mtx), rownames(mtx.BM))
# Reindex both matrices to this union, filling in zeros where necessary
mtx <- mtx[match(all_genes, rownames(mtx)), , drop=FALSE]
mtx.BM <- mtx.BM[match(all_genes, rownames(mtx.BM)), , drop=FALSE]

# Create a count matrix that only takes into account the intersecting genes
common_genes <- intersect(rownames(mtx), rownames(mtx.BM))
mtx <- mtx[match(common_genes, rownames(mtx)), , drop = FALSE]
mtx.BM <- mtx.BM[match(common_genes, rownames(mtx.BM)), , drop = FALSE]

# Replace NAs with zeros
mtx[is.na(mtx)] <- 0
mtx.BM[is.na(mtx.BM)] <- 0

#Calculate nCount_RNA and nFeature_RNA to add them to the meta_data for the custom matrix
# Filter out cells with fewer than min_gene_number genes detected
cells_to_keep <- colSums(mtx > 0) >= min_gene_number
mtx <- mtx[, cells_to_keep]

nCount_RNA <- colSums(mtx)
nFeature_RNA <- colSums(mtx > 0)

#Extract the meta data and bind it
meta.data <- l.meta.mtx[[1]]
meta.data$nCount_RNA <- nCount_RNA[rownames(meta.data)]
meta.data$nFeature_RNA <- nFeature_RNA[rownames(meta.data)]

meta.data.BM <- seurat.BM@meta.data

# The bone marrow datasets are massive compared to IRIS, so let's subset
meta.data.BM %>% count(donor_id)

selected.id <- c("TSP2", "TSP14")
meta.data.BM <- meta.data.BM %>% filter(donor_id %in% selected.id)
meta.data.BM$expID <- meta.data.BM$donor_id
meta.data.BM$cell_id <- rownames(meta.data.BM)

all_columns <- union(colnames(meta.data), colnames(meta.data.BM))
# Add missing columns to custom_meta_data with appropriate default values
for (col in setdiff(all_columns, colnames(meta.data))) {
  meta.data[[col]] <- NA
}
# Add missing columns to meta_data_10X with appropriate default values
for (col in setdiff(all_columns, colnames(meta.data.BM))) {
  meta.data.BM[[col]] <- NA
}
#Only keep those cells that passed the filtering
filtered.BM.cells <- rownames(meta.data.BM)
mtx.BM <- mtx.BM[, filtered.BM.cells] # Restrict the count matrix. Note, colnames(mtx) == rownames(meta.data)

# Note: ncol(mtx) + ncol(subsampled mtx.BM) ≈ ncol(mtx) * 1.1  because ncol(mtx.BM) 3 X ncol(mtx)
# Number of cells in reference matrix
n_ref <- ncol(mtx)
# Maximum allowed total (10% more)
target_BM_cells <- ceiling(n_ref * 1.1)

# Cap to available BM cells if needed
target_BM_cells <- min(target_BM_cells, ncol(mtx.BM))

# Subsample BM cells
set.seed(42)
subsampled_cells <- sample(colnames(mtx.BM), target_BM_cells)
mtx.BM.sub <- mtx.BM[, subsampled_cells]
meta.data.BM.sub <- meta.data.BM[subsampled_cells, ]

# Combine metadata and matrices
combined.meta.data <- rbind(meta.data, meta.data.BM)
combined.matrix <- cbind(mtx, mtx.BM)

#####
#Batch integration with Seurat
#####
#Seurat
integration.results <- seurat_integrate(
  signature_cell_cycle_human_mouse,
  species,
  c(),
  combined.matrix,
  combined.meta.data,
  integrate_over,
  min_cells_percent,
  min_gene_number,
  mito_cutoff,
  100, # dims_use
  vars_to_regress,
  split.layers,
  integrate.layers,
  integration.type,
  paste0(PATH_output_figures, "/QC_plots"))

seurat.object <- integration.results[[1]]
dims_use <- 75
integration.name <- integration.results[[3]]

resolutions <- seq(0.5, 1.5, by = 0.1) # Define the range of resolutions
find_cluster_resolution(seurat.object, 
                        dims_use, resolutions, min.pct, logfc.threshold, integration.name,
                        "Leiden", paste0(PATH_output_figures, "/UMAP_res_plots"))

if(!file.exists(paste0(PATH_output_figures, "/unnamed_plots"))){
  dir.create(paste0(PATH_output_figures, "/unnamed_plots"), recursive = TRUE)
  message(paste("Directory created: ", paste0(PATH_output_figures, "/unnamed_plots")))
} else {
  message(paste("Directory already exists: ", paste0(PATH_output_figures, "/unnamed_plots")))
}

res <- 0.9

seurat.object <- cluster_UMAP_seurat(seurat.object, res, dims_use, integration.name,
                                     reduction.name, "Leiden", 
                                     paste0(PATH_output_figures, "/unnamed_plots"))

# Save seurat.object
saveRDS(file = paste0(PATH_output_objects, "/seurat_object_integrated.Rds"),
        seurat.object)

print(paste("Fully integrated seurat object saved to: ",
            paste0(PATH_output_objects, "/seurat_object_integrated.Rds")))

# seurat.object <- readRDS(file = paste0(PATH_output_objects, "/seurat_object_integrated.Rds"))

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
  "Clustering Algorithm: ", "Leiden", "\n",
  "UMAP Resolution: ", res, "\n",
  "UMAP Integration Name: ", integration.name, "\n",
  "UMAP Reduction Name: ", reduction.name, "\n"
)

writeLines(text_content, con = paste0(PATH_output, "/integration_parameters.txt"))

#############################
#Plotting
#############################
#####
#QC
#####
generate_QC_plots(seurat.object, res, 
                  reduction.name, c("NG078", "NG075", "JP314", "TF061", "TF057", # Matched HDs
                                    "NG064", "TF051", "TF052", "JP317", "TF064", "TF066", # CHUV NPM1-mut
                                    "JP315", "TF058", "NG080", "TF063", "NG082", "NG083", "TF065", # SL NPM1-mut
                                    "JP303", "JP304", "TF050", "TF047", "TF048", # CHUV KMT2Ar
                                    "NG079", "JP318", "TF062", "TF059", "JP316", # SL KMT2Ar
                                    "JP299", "NG070"), # PB-NPM1
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
  seurat.object <- AddModuleScore(seurat.object,
                                  features = list(gene_list[[score_name]]),
                                  name = score_name)
}
for(score_name in unique(top.BMMC.DEGs$cluster)){
  print(score_name)
  gene_list[[score_name]] <- top.BMMC.DEGs[top.BMMC.DEGs$cluster == (score_name), ]$gene
  seurat.object <- AddModuleScore(seurat.object, features = list(gene_list[[score_name]]),
                                  name = score_name)
}

# Remove the "1" suffix from the new score columns in meta.data
# Renaming columns in the meta.data slot
names(seurat.object@meta.data)[(ncol(seurat.object@meta.data) - length(gene_list) + 1):ncol(seurat.object@meta.data)] <- names(gene_list)

# Remove underscores from gene names for plotting
features_clean <- gsub("_", " ", names(gene_list))  # or replace "_" with " " for spaces

mincoff <- c("q60")
plots <- FeaturePlot(seurat.object, features = names(gene_list),
                     reduction = reduction.name, cols = c("lightgrey", "deepskyblue", "blue"),
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
ggsave(filename = paste0(PATH_output_figures, "/unnamed_plots",
                         "/ct_ftr_mincoff_", mincoff, ".png"),
       plot = p, width = 32, height = 32)

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

######
#Cluster annotation
######
#Rename clusters with major cell group attribution and then subselect again

#Do not use the same names as for signature_PBMCs
new.cluster.names <- c("", "", # Cluster 1-2
                       "", "", # Cluster 3-4
                       "", "", # Cluster 5-6
                       "", "", # Cluster 7-8
                       "", "", # Cluster 9-10
                       "", "", # Cluster 11-12
                       "", "", # Cluster 13-14
                       "", "", # Cluster 15-16
                       "", "", # Cluster 17-18
                       "", "", # Cluster 19-20
                       "", "", # Cluster 21-22
                       "", "", # Cluster 23-24
                       "", "", # Cluster 25-26
                       "", "", # Cluster 27-28
                       "", "" # Cluster 29-30
                       )

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

########
# Violin plots for specific blast/oncological markers
########
#Subset the seurat object to contain only the targeted experiment
PATH_output_vln <- paste0(PATH_output_figures, "/unnamed_vln_plots/Er_prog_cells/")
dir.create(PATH_output_vln, recursive = TRUE)

#Find genes with specific sequences
# Extract all gene names from the Seurat object
all_genes <- rownames(seurat.object)
target_genes <- grep("HB", all_genes, value = TRUE, ignore.case = TRUE)

ftrs.2.plot <- c("LYL1", "DLK1", "GATA1", "HBD", "ITGA2B", "FCGR1A", "CSF1")

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

for(ftr in ftrs.2.plot){
  p <- VlnPlot(seurat.object, features = ftr) + 
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 90))
  print(p)
  ggsave(filename = paste0(PATH_output_vln, ftr, "_vln_plot.jpeg"),
         plot = p, width = 12, height = 8)
}

#########
# Volcano plot
#########
for(clst in unique(l.new.markers[[1]]$cluster)){
  mrkrs <- l.new.markers[[1]] %>% filter(cluster == clst)
  
  # The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
  mrkrs$diffexpressed <- "NO" # add a column of NAs
  
  mrkrs$diffexpressed[mrkrs$avg_log2FC > 0.6 & mrkrs$p_val_adj < 0.05] <- "UP" # If avg_log2FC > 0.6 and p_val < 0.05, set as "UP"
  mrkrs$diffexpressed[mrkrs$avg_log2FC < -0.6 & mrkrs$p_val_adj < 0.05] <- "DOWN" # If avg_log2FC < -0.6 and p_val < 0.05, set as "DOWN"
  
  mrkrs$label <- NA # Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
  mrkrs$label[mrkrs$diffexpressed != "NO"] <- mrkrs$gene[mrkrs$diffexpressed != "NO"]
  
  volcano.colors <- c("blue", "red", "black")
  names(volcano.colors) <- c("DOWN", "UP", "NO")
  
  p <- ggplot(data=mrkrs, aes(x=avg_log2FC, y=-log10(p_val_adj), 
                              col=diffexpressed, label=label)) + 
    geom_point()
  p2 <- p + geom_vline(xintercept = c(-0.6,0.6), col="red") +
    geom_hline(yintercept = -log10(0.05), col="red") + 
    scale_color_manual(values=volcano.colors) + geom_text_repel() + ggtitle(clst)
  # print(p2)
  clst_name <- gsub(" ", "_", clst)
  ggsave(filename = paste0(PATH_output_figures, "/volcano_plots/", clst_name,
                           "_volcano.png"), plot = p)
}
#########
# Gene ontology
#########
GO_path = paste0(PATH_output_figures, "/GO_plots")

# Subselect clusters if needed
filt.mrkrs.4.go <- l.new.markers[[2]] %>%
  filter(grepl("NPM", cluster, ignore.case = FALSE))

l_GO_results <- calculate_GO_terms(filt.mrkrs.4.go, species, 3, GO_path)
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

########
# Create bubble plots
########
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
cell_counts$expID <- factor(cell_counts$expID, levels = expID_order)

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
ggsave(filename = paste0(PATH_output_figures, "/named_plots/AML_bubble_plot.jpeg"), 
       plot = p, width = 12, height = 8)

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
ggsave(filename = paste0(PATH_output_figures, "/named_plots/gene_bubble_plot.jpeg"), 
       plot = p, width = 12, height = 5*ngenes)

# Plot canonical marker bubble plot per cluster
mrkrs.2.plot <- c("CD14", "KMT2A", "MLLT10", "IDH1", "NPM1", "SRSF2",
                  "PTPRC", "CD34", "KIT", "HLA.DRA", "CD38", "CD33", 
                  "CD4", "CD64", "CD300")
path <- paste0(PATH_output_figures, "/bubble_plots/KMT2Ar-a_CH5_AML_clinical_mrkrs.jpeg")

bubble.data <- l.new.markers[[1]] %>%
  filter(gene %in% mrkrs.2.plot) %>%
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
         width = 12, height = 2*length(mrkrs.2.plot))
} else {
  ggsave(filename = path, plot = p,
         width = 12, height = length(mrkrs.2.plot))
}

##########
# Create signature dataframes and lists
##########
l.new.markers <- readRDS(file = paste0(PATH_output_tables, "/l_markers.Rds"))

#Order and filter cluster DEGs
l.new.markers[[2]] %>%
  group_by(cluster) %>%
  filter(avg_log2FC >= 0.2) %>%
  top_n(n = 500, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) -> new.markers.ordered

# Subselect clusters if needed
filtered.mrkrs <- new.markers.ordered %>%
  filter(grepl("health", cluster, ignore.case = TRUE))

write.table(new.markers.ordered, file = paste0(PATH_output_tables, "/NPM1_NG064_TSP21_TSP25_sig.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
