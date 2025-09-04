rename_clusters <- function(seurat.object, new.cluster.names, reduction.name, path){
  #Set the old cluster names
  old.cluster.names <- levels(Idents(seurat.object))

  #Create a named vector for the new cluster names
  rename.vector <- setNames(new.cluster.names, old.cluster.names)

  #Rename the Idents in the seurat.object
  seurat.object <- RenameIdents(seurat.object, rename.vector)

  # Ensure the output directory exists
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  
  # Count the number of unique clusters
  cluster_labels <- unique(Idents(seurat.object))
  num_clusters <- length(cluster_labels)
  
  # Dynamically set height based on number of clusters (0.25 units per cluster)
  base_height <- 8  # Base height
  legend_height <- max(3, num_clusters * 0.25)  # Ensure a minimum extra height
  final_height <- base_height + legend_height  # Compute final height
  
  # Generate UMAP plot with vertical legend
  p <- DimPlot(seurat.object, reduction = reduction.name, label.size = 5, label = FALSE) +
    theme(
      legend.position = "bottom",  # Keep legend to the right
      legend.box = "vertical",  # Stack legend vertically
      legend.justification = "center",
      plot.margin = margin(10, 10, 10, 10)  # Adjust plot margins
    ) +
    guides(color = guide_legend(ncol = min(3, num_clusters),
                                override.aes = list(size = 4)))  # Force max 3 columns
  
  
  print(p)
  
  # Save plot with fixed width but dynamic height
  ggsave(filename = paste0(path, "/UMAP_named.png"), plot = p, 
         width = 16, height = final_height, limitsize = FALSE)
  
  
  return(seurat.object)
}

calculate_GO_terms <- function(top.DEGs, species, n_ftrs, path){
    l.top.DEGs <- list()
    for(cluster_name in unique(top.DEGs$cluster)){
    l.top.DEGs[[cluster_name]] <- top.DEGs %>% filter(cluster == cluster_name)
    }
    l_GO_results <- gene_ontology_compute(species, l.top.DEGs)

    l_top5_GO <- lapply(l_GO_results, function(go_result) {
    go_result %>%
        arrange(`Rank in weight`) %>% # Adjust column name to your p-value adjustment method
        head(n_ftrs)
    })

    # Ensure the output directory exists
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
    }

    p <- gene_ontology_heatmap(l_top5_GO, 0.025, n_ftrs)
    print(p)
    ggsave(filename = paste0(path, "/GO_heatmap.png"),
       plot = p$ggplot_heat_GO, width = 8, height = 12)

    return(l_GO_results)
}

# Function to order and plot GO terms for each cluster
plot_go_terms_for_cluster <- function(cluster_results, cluster_name, 
                                      top_n_terms, path) {
    # Convert weight column to numeric
    cluster_results$weight <- as.numeric(as.character(cluster_results$weight))
    
    # Filter for significant terms based on the weight p-value
    significant_terms <- cluster_results %>%
        filter(weight < 0.05) %>%
        arrange(weight)
    
    # Select the top N most significant terms for plotting
    top_terms <- significant_terms %>%
        top_n(-top_n_terms, weight)
    
    # Reorder terms for plotting
    top_terms <- top_terms %>%
      mutate(Term = str_wrap(Term, width = 50)) %>% # Wrap after 50 characters
      mutate(Term = factor(Term, levels = unique(Term[order(weight)])))
    
    # Wrap terms so that GO terms are printed in entirety
    top_terms <- top_terms 
    
    # Plotting the top terms
    plot <- ggplot(top_terms, aes(x = Term, y = -log10(weight))) +
        geom_bar(stat = "identity", fill = "skyblue") +
        coord_flip() +
        labs(title = paste("Top GO Terms for Cluster", cluster_name),
            x = "GO Term",
            y = "-log10(Weight P-Value)") +
        theme_bw() +  # Use theme_bw for white background
        theme(
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white", color = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8), # Adjust for size readability
        plot.margin = margin(t = 10, r = 30, b = 10, l = 10) # Add right side padding
        )
    
    # Save the plot
    ggsave(paste0(path, "/GO_", cluster_name, ".png"), plot, width = 10, height = 8)
    
    return(plot)
}


#Select GO terms for pathways linked to specific terms
#Function to filter GO terms based on keywords
filter_GO_terms <- function(l_GO_results, keywords, n_ftrs) {
  # Use lapply to process each element (GO results for each cluster)
  l_filtered_GO <- lapply(l_GO_results, function(go_result) {
    #Check that the 'Term' column exists
    if (!"Term" %in% names(go_result)) {
      warning("Column 'Term' not found in one element of l_GO_results")
      return(NULL)
    }
    
    #Initialize an empty dataframe for concatenating results
    filtered_results <- data.frame()
    
    #Loop through each keyword and subset the data
    for(keyword in keywords) {
      keyword_subset <- go_result %>%
        filter(grepl(keyword, Term, ignore.case = TRUE)) %>%
        mutate(Keyword_Matched = keyword) #Track the matched keyword
      
      #Concatenate the subsetted results
      filtered_results <- bind_rows(filtered_results, keyword_subset)
    }
    
    #Remove duplicate rows (in case multiple keywords match the same terms)
    filtered_results <- filtered_results %>%
      distinct(GO.ID, Term, .keep_all = TRUE) # Keep only one instance of each GO.ID and Term
    
    #Order the concatenated results by descending `Rank in weight`
    filtered_results <- filtered_results %>%
      arrange(`Rank in weight`) # Sorting in descending order
    
    return(filtered_results)
  })
  
  l_top5_GO <- lapply(l_filtered_GO, function(go_result) {
    go_result %>%
      arrange(`Rank in weight`) %>% # Adjust column name to your p-value adjustment method
      head(n_ftrs)
  })
  
  # Ensure the output directory exists
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  
  p <- gene_ontology_heatmap(l_top5_GO, 0.025, n_ftrs)
  print(p)
  #ggsave(filename = paste0(path, "/GO_heatmap.png"),
  #       plot = p$ggplot_heat_GO, width = 8, height = 12)
  
  return(list(l_filtered_GO, l_top5_GO))
}