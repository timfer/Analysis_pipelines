annotate_doublets <- function(PATH_input_IRIS_imaging, sample_ids_IRIS){
  l_image_doublet_ROI_annotation <- list()
  for(i in 1:length(sample_ids_IRIS)){
    sample_ids_IRIS_i <- sample_ids_IRIS[i]
    PATH_input_IRIS_imaging_ROIcount <- paste0(PATH_input_IRIS_imaging,
                                               "/",sample_ids_IRIS_i,
                                               "/metadata/cell_count_metadata.csv")
                                               # "/image_annotation/algo/cells_per_ROI/complete_df.csv")
    t_image_doublet_ROI_annotation_i <- read.csv(PATH_input_IRIS_imaging_ROIcount)
    l_image_doublet_ROI_annotation[[i]] <- t_image_doublet_ROI_annotation_i
  }
  names(l_image_doublet_ROI_annotation) <- sample_ids_IRIS
  lapply(l_image_doublet_ROI_annotation, nrow)
  return(l_image_doublet_ROI_annotation)
}

create_matrix_meta_data <- function(sample_ids_IRIS, l.image.doublet.ROI.annotation){
  l_pivot_table <- list()
  l_UMI_all_cells_SYMB <- list()
  l_SM_annotation <- list()
  pivot_full_id <- "Full_Cell_ID"
  meta_full_id <- "full_cell_id"
  for(i in 1:length(sample_ids_IRIS)){
    print(i)
    print(sample_ids_IRIS[i])
    
    #Load files
    load(paste0(PATH_input_IRIS_sequencing,"/",sample_ids_IRIS[i],"/pipeline_output/",
                species,"/analysis/objects/objects_plotting.RData"), 
         envir = parent.frame(), verbose = FALSE)
    #Select species mixing annotation file
    SM_file_path <- paste0(PATH_input_IRIS_sequencing,"/",
                           sample_ids_IRIS[i],"/results/tables/")
    pattern <- paste0(sample_ids_IRIS[i], "_transcriptome_species_annotation.*\\.txt$")
    SM_file <- list.files(path = SM_file_path, pattern = pattern, full.names = TRUE)
    
    #Select cells of species choice
    SM_annotation_i <- read.delim(SM_file)
    # Find the column regardless of case
    # pivot_full_id <- names(SM_annotation_i)[match(tolower("full_cell_id"), tolower(names(SM_annotation_i)))]
    
    Cells_wanted_i <- subset(SM_annotation_i, species_annotation == species)[[pivot_full_id]]
    print(paste0(nrow(pivot_table), "| cells detected"))
    
    #Pivot table
    pivot_table_i <- subset(pivot_table, pivot_table[[pivot_full_id]] %in% Cells_wanted_i)
    # This pivot table tracks human cells but not single cells
    print(paste0(nrow(pivot_table_i), "| human cells"))
    
    #Select cells with min_UMI
    pivot_table_i <- subset(pivot_table_i, n_UMI > min_nUMI)
    print(paste0(nrow(pivot_table_i), "| cells >3000UMI"))
    
    #Select encapsulation with only one cell in ROI
    t_cells_per_ROI_i <- l.image.doublet.ROI.annotation[[i]] # This metadata table only tracks single cells
    # meta_full_id <- names(t_cells_per_ROI_i)[match(tolower("full_cell_id"), tolower(names(t_cells_per_ROI_i)))]
    
    t_cells_per_ROI_i <- subset(t_cells_per_ROI_i, t_cells_per_ROI_i[[meta_full_id]] %in% 
                                  pivot_table_i[[pivot_full_id]])
    Cells_wanted_i <- subset(t_cells_per_ROI_i, cell_count == "1")[[meta_full_id]]
    
    #Pivot table
    pivot_table_i <- subset(pivot_table_i, pivot_table_i[[pivot_full_id]] %in% Cells_wanted_i)
    print(paste0(nrow(pivot_table_i), "| cells single encapsulations"))
    
    #Count matrix
    ncol(UMI_all_cells_SYMB)
    UMI_all_cells_SYMB_i <- UMI_all_cells_SYMB[,colnames(UMI_all_cells_SYMB) %in% 
                                                 pivot_table_i[[pivot_full_id]]]
    ncol(UMI_all_cells_SYMB_i)                     
    #Store data
    l_SM_annotation[[i]] <- SM_annotation_i
    l_pivot_table[[i]] <- pivot_table_i
    l_UMI_all_cells_SYMB[[i]] <- UMI_all_cells_SYMB_i
  }
  #names to lists
  names(l_SM_annotation) <- sample_ids_IRIS
  names(l_pivot_table) <- sample_ids_IRIS
  names(l_UMI_all_cells_SYMB) <- sample_ids_IRIS
  lapply(l_pivot_table, nrow)
  lapply(l_pivot_table, ncol)
  lapply(l_UMI_all_cells_SYMB, nrow)
  lapply(l_UMI_all_cells_SYMB, ncol)
  pivot_all <- do.call(rbind, l_pivot_table)
  
  #QC on a glance
  print("meadian(mappedRead)")
  print(lapply(l_pivot_table, function(x){median(x$n_mapped_reads)}))
  print("meadian(nUMI)")
  print(lapply(l_pivot_table, function(x){median(x$n_UMI)}))
  print("mean(nUMI)")
  print(lapply(l_pivot_table, function(x){mean(x$n_UMI)}))
  print("sd(nUMI)")
  print(lapply(l_pivot_table, function(x){sd(x$n_UMI)}))
  print("median(nGenes)")
  print(lapply(l_pivot_table, function(x){median(x$n_genes)}))
  print("nCells")
  print(unlist(lapply(l_pivot_table, nrow)))
  
  #####
  #Prep meta data
  #####
  #Create meta data
  lapply(l_pivot_table, nrow)
  meta_data <- do.call(rbind, l_pivot_table)
  rownames(meta_data) <- meta_data[[pivot_full_id]]
  names(meta_data)[names(meta_data) == "Sample_Plate"] <- "expID"
  levels(as.factor(meta_data$expID))

  #####
  #Prep matrix
  #####
  #https://stackoverflow.com/questions/7739578/merge-data-frames-based-on-rownames-in-r
  #test <- merge(l_UMI_all_cells_SYMB[[1]], l_UMI_all_cells_SYMB[[2]], by=0, all=TRUE)
  #Merge matrices before integrating data, this should limit integration artifacts
  matrix <- l_UMI_all_cells_SYMB[[1]]
  if(length(l_UMI_all_cells_SYMB)>1){
    for(i in 1:(length(l_UMI_all_cells_SYMB)-1)){
      print(names(l_UMI_all_cells_SYMB)[i])
      print(names(l_UMI_all_cells_SYMB)[i+1])
      
      exp_name <- names(l_UMI_all_cells_SYMB)[i]
      exp_name_2 <- names(l_UMI_all_cells_SYMB)[i+1]
      
      remaining.IDs <- l_pivot_table[[exp_name]][[pivot_full_id]]
      remaining.IDs.2 <- l_pivot_table[[exp_name_2]][[pivot_full_id]]
      
      l_UMI_all_cells_SYMB[[i]] <- l_UMI_all_cells_SYMB[[i]][, colnames(l_UMI_all_cells_SYMB[[i]]) %in%
                                                              remaining.IDs]
      l_UMI_all_cells_SYMB[[i+1]] <- l_UMI_all_cells_SYMB[[i+1]][, colnames(l_UMI_all_cells_SYMB[[i+1]]) %in%
                                                                  remaining.IDs.2]
      
      matrix <- merge(matrix, l_UMI_all_cells_SYMB[[i+1]], by=0, all=TRUE)
      rownames(matrix) <- matrix$Row.names
      matrix <- matrix[,2:ncol(matrix)]
      rownames(matrix)[grep("^MT.", rownames(matrix))]
    }
  }
  sum(unlist(lapply(l_UMI_all_cells_SYMB, ncol)))
  ncol(matrix)
  nrow(matrix)
  #Fill NAs
  matrix[is.na(matrix)] <- 0 #Set all NAs in the matrix to 0
  matrix <- as.matrix(matrix)

  #Fix "MT-": Reassign MT names to mitochondrial genes
  rownames(matrix)[grep("^MT\\.", rownames(matrix))]
  rownames(matrix) <- gsub("^MT\\.", "MT-", rownames(matrix))
  rownames(matrix)[grep("^MT-", rownames(matrix))]

  return(list(meta_data, matrix))
}