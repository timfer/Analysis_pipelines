#Scope | Support functions for IRIS scRNA-seq, imaging and seurat formats
# 1) Combine seurat-meta data and IRIS scRNA-seq meta-data from pivot_table
# 2) Combine seurat-meta data with IRIS imaging meta-data
# 3) GO analysis on DEGs genes
#Script-Type | Function
#Input | Various (details see per function)
#Output | Various (details see per function)
#Author | Joern Pezoldt
#Date | init 30.12.2023




#####
# 1) Combine seurat-meta data and IRIS scRNA-seq meta-data from pivot_table
#####
#RETURNs: seurat.object with replaced @meta.data containing data from original seurat.object and the pivot-table

add_pivot_to_seurat <- function(
    # Seurat object from  integrated seurat object via seurat_batch()
    seurat.object,
    # pivot table of all integrated datasets via rbind()
    pivot_table_all){
  #Store meta data
  MetaData <- seurat.object@meta.data
  order_meta <- rownames(MetaData)
  IDs_LiD_in_Seurat <- rownames(MetaData)
  #Make compatible table order
  MetaData <- MetaData[order(rownames(MetaData)),]
  MetaData$Full_Cell_ID <- rownames(MetaData)
  pivot_all_in_Seurat <- subset(pivot_table_all, Full_Cell_ID %in% IDs_LiD_in_Seurat)
  pivot_all_in_Seurat <- pivot_all_in_Seurat[order(pivot_all_in_Seurat$Full_Cell_ID),]
  #merge
  meta_seurat_iris <- merge(MetaData, pivot_all_in_Seurat, by = "Full_Cell_ID",
                            all.x = TRUE)
  meta_seurat_iris <- meta_seurat_iris[match(order_meta, meta_seurat_iris$Full_Cell_ID),]
  #Store IDs
  Full_Cell_ID <- meta_seurat_iris$Full_Cell_ID
  meta_seurat_iris$Full_Cell_ID <- NULL
  #Add IDs at relevant sections in table
  rownames(meta_seurat_iris) <- Full_Cell_ID
  meta_seurat_iris$Full_Cell_ID <- Full_Cell_ID
  #Save to seurat objects
  seurat.object@meta.data <- meta_seurat_iris
  #return
  seurat.object
}


#####
# 2.1) Combine seurat-meta data with cell count ROI
#####
#RETURNs: seurat.object with replaced @meta.data containing data from original seurat.object and the numbers of cells in imaging area

add_ROIcellCOUNT_to_seurat <- function(
    # Seurat object from  integrated seurat object via seurat_batch()
  seurat.object,
  # cell count in imaging area/ROI per dataset via rbind()
  t_count_cells_ROI){
  #Store meta data
  MetaData <- seurat.object@meta.data
  order_meta <- rownames(MetaData)
  IDs_LiD_in_Seurat <- rownames(MetaData)
  #Make compatible table order
  MetaData <- MetaData[order(rownames(MetaData)),]
  MetaData$Full_Cell_ID <- rownames(MetaData)
  cellcounts_in_Seurat <- subset(t_count_cells_ROI, Full_Cell_ID %in% IDs_LiD_in_Seurat)
  cellcounts_in_Seurat <- cellcounts_in_Seurat[order(cellcounts_in_Seurat$Full_Cell_ID),]
  cellcounts_in_Seurat <- cellcounts_in_Seurat[,c("Full_Cell_ID","cell_count")]
  #merge
  meta_seurat_iris <- merge(MetaData, cellcounts_in_Seurat, by = "Full_Cell_ID",
                            all.x = TRUE)
  meta_seurat_iris <- meta_seurat_iris[match(order_meta, meta_seurat_iris$Full_Cell_ID),]
  #Store IDs
  Full_Cell_ID <- meta_seurat_iris$Full_Cell_ID
  meta_seurat_iris$Full_Cell_ID <- NULL
  #Add IDs at relevant sections in table
  rownames(meta_seurat_iris) <- Full_Cell_ID
  meta_seurat_iris$Full_Cell_ID <- Full_Cell_ID
  v_cell_count <- meta_seurat_iris$cell_count
  v_cell_count[is.na(v_cell_count)] <- "no_image"
  meta_seurat_iris$cell_count <- v_cell_count
  #Save to seurat objects
  seurat.object@meta.data <- meta_seurat_iris
  #return
  seurat.object
}

#####
# 2.2) Combine seurat-meta data with features per cell
#####
#RETURNs: seurat.object with replaced @meta.data containing data from original seurat.object and features computed per cell
add_image_features_to_seurat <- function(
    # Seurat object from  integrated seurat object via seurat_batch()
  seurat.object,
  # cell count in imaging area/ROI per dataset via rbind()
  t_imaging_features){
  #Store meta data
  MetaData <- seurat.object@meta.data
  order_meta <- rownames(MetaData)
  IDs_LiD_in_Seurat <- rownames(MetaData)
  #Make compatible table order
  MetaData <- MetaData[order(rownames(MetaData)),]
  MetaData$Full_Cell_ID <- rownames(MetaData)
  features_in_Seurat <- subset(t_imaging_features, Full_Cell_ID %in% IDs_LiD_in_Seurat)
  features_in_Seurat <- features_in_Seurat[order(features_in_Seurat$Full_Cell_ID),]
  #merge
  meta_seurat_iris <- merge(MetaData, features_in_Seurat, by = "Full_Cell_ID",
                            all.x = TRUE)
  meta_seurat_iris <- meta_seurat_iris[match(order_meta, meta_seurat_iris$Full_Cell_ID),]
  #Store IDs
  Full_Cell_ID <- meta_seurat_iris$Full_Cell_ID
  meta_seurat_iris$Full_Cell_ID <- NULL
  #Add IDs at relevant sections in table
  rownames(meta_seurat_iris) <- Full_Cell_ID
  meta_seurat_iris$Full_Cell_ID <- Full_Cell_ID
  meta_seurat_iris[is.na(meta_seurat_iris)] <- "not_in_ROI"
  #Save to seurat objects
  seurat.object@meta.data <- meta_seurat_iris
  #return
  seurat.object
}





##########################################
#Gene Ontology
##########################################
#####
# 3.1) GO analysis on DEGs genes
#####
#RETRURNs: list of GOs per list of differentially expressed genes suplied with p-values, GO:IDs, names...

gene_ontology_compute <- function(
  #species
  species,
  #named list of differentially expressed genes
  ## each list element contains one table with at least columns "gene", "p_val"
  list_DEGs){
      #Load libraries
      library(stringr)
      library(ggplot2)
      library(foreign)
      library(GO.db)
      library(topGO)
      library(org.Hs.eg.db)
      library(org.Mm.eg.db)

      #####
      #Make a gene2GO list
      #####
      if(species == "human"){
        x <- org.Hs.egGO}
      if(species == "mouse"){
        x <- org.Mm.egGO}
      # Get the entrez gene identifiers that are mapped to a GO ID
      mapped_genes <- mappedkeys(x)
      # Build a list GeneID and GOs
      GeneID2GO <- list()
  
      xx <- as.list(x[mapped_genes])
      
      for(i in 1:length(xx)){
        #Initiate vector to collect GOIDs for geneID i
        GO_vector <- c()
        #Get geneID
        Gene_ID <- ls(xx[i])
        #grab data for geneID_i
        temp <- xx[[i]]
        
        #check Ontology category
        for(i in 1:length(temp)){
          category <- as.character(temp[[i]]["Ontology"])
          #if ontology category matches collect GOIDs  (flex)
          if(category == "BP"){
            temp_GOID <- as.character(temp[[i]]["GOID"])
            GO_vector <- c(GO_vector, temp_GOID)
          }
          #print(GO_vector)
        }
        #Generate list name geneID, content
        GO_IDs <- GO_vector 
        GeneID2GO[[Gene_ID]] <- GO_IDs 
      }
  
      GO2GeneID <- inverseList(GeneID2GO)

      #####
      #Setup input for GO analysis
      #####
      l_input_DEGs <- list_DEGs
      new.cluster.ids <- names(l_input_DEGs)
      myfiles <- l_input_DEGs
      #GeneSymbol to GeneID
      
      if(species == "human"){
        for(i in 1:length(myfiles)){
          #i =1
          idfound <- myfiles[[i]]$gene %in% mappedRkeys(org.Hs.egSYMBOL)
          SYMBOL <- toTable(org.Hs.egSYMBOL)
          head(SYMBOL)
          m <- match(myfiles[[i]]$gene, SYMBOL$symbol)
          GENE_ID <- SYMBOL$gene_id[m]
          store_i <- cbind(GENE_ID, myfiles[[i]])
          colnames(store_i) <- c("GENE_ID",colnames(store_i)[2:ncol(store_i)])
          #store_i <- store_i[complete.cases(store_i), ]
          myfiles[[i]] <- store_i
        }
      }
      if(species == "mouse"){
        for(i in 1:length(myfiles)){
          #i =1
          idfound <- myfiles[[i]]$gene %in% mappedRkeys(org.Mm.egSYMBOL)
          SYMBOL <- toTable(org.Mm.egSYMBOL)
          head(SYMBOL)
          m <- match(myfiles[[i]]$gene, SYMBOL$symbol)
          GENE_ID <- SYMBOL$gene_id[m]
          store_i <- cbind(GENE_ID, myfiles[[i]])
          colnames(store_i) <- c("GENE_ID",colnames(store_i)[2:ncol(store_i)])
          #store_i <- store_i[complete.cases(store_i), ]
          myfiles[[i]] <- store_i
        }
      }
      
      #Perform GO analysis single
      #determine number of clusters: split by cluster
      #link gene IDs with p-values
      #perform GO
      l_gene_pval <- list()
      l_all_gene_pval <- list()
      
      for(i in 1:length(myfiles)){
        #i= 1
        genes_of_interest_i <- as.character(myfiles[[i]]$GENE_ID)
        genes_pval_i <- myfiles[[i]]$p_val
        names(genes_pval_i) <- genes_of_interest_i
        l_gene_pval[[i]] <- genes_pval_i 
      }
      print(class(l_gene_pval))
      print(length(l_gene_pval))
      l_all_gene_pval[[length(l_all_gene_pval)+1]] <- l_gene_pval
      
      
      #get one gene list
      l_gene_pval <- l_all_gene_pval[[1]]
      
      #Gene universe
      #Genes in DEGs list
      geneNames <- unique(unlist(lapply(myfiles, function(x){x$GENE_ID})))
      #Size gene universe
      length(geneNames)

      #gene lists
      l_gene_List <- list()
      for(i in 1:length(l_gene_pval)){
        geneList_i <- factor(as.integer(geneNames %in% names(l_gene_pval[[i]])))
        names(geneList_i) <- geneNames
        l_gene_List[[i]] <- geneList_i 
      }
      length(l_gene_List)
      
      #####
      #GO analysis
      #####
      List_allRes <- list()
      #Do  GO statistics for all gene lists
      for(i in 1:length(l_gene_List)){
        #Access the gene lists and p-values for the differentially expressed genes
        #i = 1
        print(paste0(i, " ||| ", names(l_input_DEGs)[i]))
        geneList <- l_gene_List[[i]]
        pvalue_of_interest <- l_gene_pval[[i]]
        
        #build GOdata object
        GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel = pvalue_of_interest,
                      annot = annFUN.gene2GO , gene2GO = GeneID2GO, nodeSize = 5)
        
        #get number of differentially expressed genes in the GOdata object
        sg <- sigGenes(GOdata)
        numSigGenes(GOdata)
        #get the number of GO_IDs that are within the applied GeneUniverse
        graph(GOdata)
        number_GOIDs <- usedGO(GOdata)
        number_nodes <- length(number_GOIDs)
        
        
        #Run statistics
        #Fisher's works with weight
        test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
        resultWeight <- getSigGroups(GOdata, test.stat)
        #KS works with elim but not with weight
        #test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
        #resultElim <- getSigGroups(GOdata, test.stat)
        #runTest a high level interface for testing Fisher
        resultFis <- runTest(GOdata, statistic = "fisher")
        #Kolmogorov-Smirnov includes p-values
        #test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
        #resultKS <- getSigGroups(GOdata,  test.stat)
        #runTest a high level interface for testing KS
        #elim.ks <- runTest(GOdata, algorithm = "elim", statistic = "ks")
        
        #make table
        allRes <- GenTable(GOdata, classic = resultFis, weight = resultWeight,
                           orderBy = "weight", topNodes = number_nodes)
        #KS = resultKS, 
        #ranksOf = "KS", 
        #make list of result tables
        List_allRes[[i]] <- allRes
      }
      names(List_allRes) <- new.cluster.ids
      #return
      List_allRes
}

#####
# 3.2) GO analysis on DEGs genes with all detected genes as universe
#####
#RETRURNs: list of GOs per list of differentially expressed genes suplied with p-values, GO:IDs, names...

gene_ontology_universe <- function(
    #species
  species,
  gene_universe,
  #named list of differentially expressed genes
  ## each list element contains one table with at least columns "gene", "p_val"
  list_DEGs){
  #Load libraries
  library(stringr)
  library(ggplot2)
  library(foreign)
  library(GO.db)
  library(topGO)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  
  #####
  #Make a gene2GO list
  #####
  if(species == "human"){
    x <- org.Hs.egGO}
  if(species == "mouse"){
    x <- org.Mm.egGO}
  # Get the entrez gene identifiers that are mapped to a GO ID
  mapped_genes <- mappedkeys(x)
  # Build a list GeneID and GOs
  GeneID2GO <- list()
  
  xx <- as.list(x[mapped_genes])
  
  for(i in 1:length(xx)){
    #Initiate vector to collect GOIDs for geneID i
    GO_vector <- c()
    #Get geneID
    Gene_ID <- ls(xx[i])
    #grab data for geneID_i
    temp <- xx[[i]]
    
    #check Ontology category
    for(i in 1:length(temp)){
      category <- as.character(temp[[i]]["Ontology"])
      #if ontology category matches collect GOIDs  (flex)
      if(category == "BP"){
        temp_GOID <- as.character(temp[[i]]["GOID"])
        GO_vector <- c(GO_vector, temp_GOID)
      }
      #print(GO_vector)
    }
    #Generate list name geneID, content
    GO_IDs <- GO_vector 
    GeneID2GO[[Gene_ID]] <- GO_IDs 
  }
  
  GO2GeneID <- inverseList(GeneID2GO)
  
  #####
  #Setup input for GO analysis
  #####
  l_input_DEGs <- list_DEGs
  new.cluster.ids <- names(l_input_DEGs)
  myfiles <- l_input_DEGs
  mygenes <- gene_universe
  #GeneSymbol to GeneID
  
  if(species == "human"){
    for(i in 1:length(myfiles)){
      #i =1
      #Translate DEGs
      idfound <- myfiles[[i]]$gene %in% mappedRkeys(org.Hs.egSYMBOL)
      SYMBOL <- toTable(org.Hs.egSYMBOL)
      head(SYMBOL)
      m <- match(myfiles[[i]]$gene, SYMBOL$symbol)
      GENE_ID <- SYMBOL$gene_id[m]
      store_i <- cbind(GENE_ID, myfiles[[i]])
      colnames(store_i) <- c("GENE_ID",colnames(store_i)[2:ncol(store_i)])
      #store_i <- store_i[complete.cases(store_i), ]
      myfiles[[i]] <- store_i
    }
    #Translate detected genes
    idfound_mygenes <- mygenes %in% mappedRkeys(org.Hs.egSYMBOL)
    SYMBOL <- toTable(org.Hs.egSYMBOL)
    head(SYMBOL)
    m <- match(mygenes, SYMBOL$symbol)
    GENE_ID <- SYMBOL$gene_id[m]
    store <- cbind(GENE_ID, mygenes)
    colnames(store) <- c("GENE_ID",colnames(store)[2:ncol(store)])
    mygenes <- store
  }
  if(species == "mouse"){
    for(i in 1:length(myfiles)){
      #i =1
      idfound <- myfiles[[i]]$gene %in% mappedRkeys(org.Mm.egSYMBOL)
      SYMBOL <- toTable(org.Mm.egSYMBOL)
      head(SYMBOL)
      m <- match(myfiles[[i]]$gene, SYMBOL$symbol)
      GENE_ID <- SYMBOL$gene_id[m]
      store_i <- cbind(GENE_ID, myfiles[[i]])
      colnames(store_i) <- c("GENE_ID",colnames(store_i)[2:ncol(store_i)])
      #store_i <- store_i[complete.cases(store_i), ]
      myfiles[[i]] <- store_i
    }
    #Translate detected genes
    idfound_mygenes <- mygenes %in% mappedRkeys(org.Mm.egSYMBOL)
    SYMBOL <- toTable(org.Mm.egSYMBOL)
    head(SYMBOL)
    m <- match(mygenes, SYMBOL$symbol)
    GENE_ID <- SYMBOL$gene_id[m]
    store <- cbind(GENE_ID, mygenes)
    colnames(store) <- c("GENE_ID",colnames(store)[2:ncol(store)])
    mygenes <- store
  }
  
  #Perform GO analysis single
  #determine number of clusters: split by cluster
  #link gene IDs with p-values
  #perform GO
  l_gene_pval <- list()
  l_all_gene_pval <- list()
  
  for(i in 1:length(myfiles)){
    #i= 1
    genes_of_interest_i <- as.character(myfiles[[i]]$GENE_ID)
    genes_pval_i <- myfiles[[i]]$p_val
    names(genes_pval_i) <- genes_of_interest_i
    l_gene_pval[[i]] <- genes_pval_i 
  }
  print(class(l_gene_pval))
  print(length(l_gene_pval))
  l_all_gene_pval[[length(l_all_gene_pval)+1]] <- l_gene_pval
  
  
  #get one gene list
  l_gene_pval <- l_all_gene_pval[[1]]
  
  #Gene universe
  #Genes in DEGs list
  geneNames <- unique(as.data.frame(mygenes)$GENE_ID)[2:length(mygenes)]
  #Size gene universe
  length(geneNames)
  
  #gene lists
  l_gene_List <- list()
  for(i in 1:length(l_gene_pval)){
    geneList_i <- factor(as.integer(geneNames %in% names(l_gene_pval[[i]])))
    names(geneList_i) <- geneNames
    l_gene_List[[i]] <- geneList_i 
  }
  length(l_gene_List)
  
  #####
  #GO analysis
  #####
  List_allRes <- list()
  #Do  GO statistics for all gene lists
  for(i in 1:length(l_gene_List)){
    #Access the gene lists and p-values for the differentially expressed genes
    #i = 1
    print(paste0(i, " ||| ", names(l_input_DEGs)[i]))
    geneList <- l_gene_List[[i]]
    pvalue_of_interest <- l_gene_pval[[i]]
    
    #build GOdata object
    GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel = pvalue_of_interest,
                  annot = annFUN.gene2GO , gene2GO = GeneID2GO, nodeSize = 5)
    
    #get number of differentially expressed genes in the GOdata object
    sg <- sigGenes(GOdata)
    numSigGenes(GOdata)
    #get the number of GO_IDs that are within the applied GeneUniverse
    graph(GOdata)
    number_GOIDs <- usedGO(GOdata)
    number_nodes <- length(number_GOIDs)
    
    
    #Run statistics
    #Fisher's works with weight
    test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
    resultWeight <- getSigGroups(GOdata, test.stat)
    #KS works with elim but not with weight
    #test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
    #resultElim <- getSigGroups(GOdata, test.stat)
    #runTest a high level interface for testing Fisher
    resultFis <- runTest(GOdata, statistic = "fisher")
    #Kolmogorov-Smirnov includes p-values
    #test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
    #resultKS <- getSigGroups(GOdata,  test.stat)
    #runTest a high level interface for testing KS
    #elim.ks <- runTest(GOdata, algorithm = "elim", statistic = "ks")
    
    #make table
    allRes <- GenTable(GOdata, classic = resultFis, weight = resultWeight,
                       orderBy = "weight", topNodes = number_nodes)
    #KS = resultKS, 
    #ranksOf = "KS", 
    #make list of result tables
    List_allRes[[i]] <- allRes
  }
  names(List_allRes) <- new.cluster.ids
  #return
  List_allRes
}


#####
#Plot Heatmaps of GO statistics
#####
#RETURNs: 
## 1) heatmap
## 2) table with combined statistics per GO ID and DEG list identifer
gene_ontology_heatmap <- function(
  #output gene_ontology_compute
  list_GO_per_cluster,
  #p-value cut-off to include in heatmap
  pval_cutoff,
  #Cutoff to plot on heatmap to condense color space
  max_pval){
  #make RETURN list
  list_return <- list()
  #load library
  library(pheatmap)
  #compile individual tables
  List_TopGOs <- list()
  for(i in 1:length(list_GO_per_cluster)){
    #i = 1
    table_i <- list_GO_per_cluster[[i]]
    table_i$weight <- as.numeric(table_i$weight)
    table_tophits <- subset(table_i, weight < pval_cutoff)
    #naive 5.0e-03
    #print(i)
    #print(nrow(table_tophits))
    List_TopGOs[[i]] <- table_tophits
  }
  #collect all TopGos
  TopGOs_vector <- c()
  for(i in 1:length(List_TopGOs)){
    table_i <- List_TopGOs[[i]]
    TopGOs_i <- as.character(table_i$GO.ID)
    TopGOs_vector <- c(TopGOs_vector, TopGOs_i) 
  }
  
  #Condense tables and sort by GO.ID
  l_topGO <- list()
  for(i in 1:length(list_GO_per_cluster)){
    #i = 1
    table_i <- list_GO_per_cluster[[i]]
    head(table_i)
    table_i_subset <- subset(table_i, table_i$GO.ID %in% TopGOs_vector)
    table_i_subset_GO_weight <- table_i_subset[,c("GO.ID","weight", "Term")]
    cluster_weight_name <- names(list_GO_per_cluster)[i]
    colnames(table_i_subset_GO_weight)[2] <- cluster_weight_name
    table_i_subset_GO_weight[,2] <- as.numeric(table_i_subset_GO_weight[,2])
    table_i_subset_GO_weight <- table_i_subset_GO_weight[order(table_i_subset_GO_weight$"GO.ID", decreasing=TRUE), ]
    l_topGO[[i]] <- table_i_subset_GO_weight
  }
  #merger
  data_TopGO_weight = Reduce(function(...) merge(..., all=T), l_topGO)
  
  ####
  #Make GO comparison heatmap
  ####
  data_heatmap <- data_TopGO_weight
  rownames(data_heatmap) <- paste0(data_heatmap$GO.ID, "_", data_heatmap$Term)
  data_heatmap[is.na(data_heatmap)] <- 1
  data_heatmap <- data_heatmap[,3:ncol(data_heatmap)]
  #log transfrom
  data_heatmap <- -log10(data_heatmap)
  #Common Minimum at 5
  data_heatmap[data_heatmap > max_pval] <- max_pval
  data_heatmap_matrix <- data.matrix(data_heatmap)
  #colnames(data_heatmap_matrix) <- new.cluster.ids
  title <- c("Differential GO scRNASeq")
  ggplot_GO_heatmap <- pheatmap(data_heatmap_matrix, cluster_rows = TRUE, legend = TRUE,
           treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = FALSE,
           scale = "none", border_color = "black", cellwidth = 10,
           cellheigth = 10, color = colorRampPalette(c("white", "grey","deepskyblue2"), space="rgb")(128),
           main = title)
  #return
  list_return[[1]] <- data_TopGO_weight
  list_return[[2]] <- ggplot_GO_heatmap
  names(list_return) <- c("table_GO_weigths","ggplot_heat_GO")
  list_return
}
