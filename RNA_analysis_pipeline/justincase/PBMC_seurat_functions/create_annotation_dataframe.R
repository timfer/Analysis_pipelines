
create_annotation_dataframe <- function(seurat.object, l.image.doublet.ROI.annotation, path){
  meta_data <- seurat.object@meta.data #Select seurat object meta.data
  meta_data$cell_type <- Idents(seurat.object) #Add cluster names as a column
  # image.annotation <- bind_rows(l.image.doublet.ROI.annotation) #Create dataframe for centroid, channel and cell_count
  # 
  # "Correct the joining"
  # merged_table <- image.annotation %>%
  #   inner_join(meta_data %>% dplyr::select(Full_Cell_ID, Cluster, 
  #                                          nCount_RNA, nFeature_RNA, n_UMI, 
  #                                          n_genes, n_all_reads, n_mapped_reads,
  #                                          percent.mt, percent.Rpl, percent.Hsp), 
  #              by = "Full_Cell_ID")
  # merged_table$nFeature_RNA <- merged_table$nFeature_RNA.y
  # merged_table$nCount_RNA <- merged_table$nCount_RNA.y
  # merged_table <- merged_table %>%
  #   mutate(Cell_subclass = Cluster)
  # merged_table <- merged_table %>%
  #   dplyr::select(-Cluster, -nFeature_RNA.x, -nFeature_RNA.y, -nCount_RNA.x, -nCount_RNA.y)
  # 
  # merged_table$Cell_subclass <- gsub(" ", "_", merged_table$Cell_subclass)
  # 
  # write.csv(merged_table, file = path, row.names = FALSE)
  write.csv(meta_data, file = path, row.names = FALSE)
}