#' Get Final Clusters
#'
#' This function updates the Idents of the clusters in the Seurat object based on the object names ending with _Clean or _Extra.
#'
#' @param SeuData Seurat object containing the data
#' @return Updated Seurat object with the Idents and celltype columns updated
#'
#' @export
#'
#' @import Seurat
#' @importFrom SeuratObject Idents
#' @importFrom stats cov
#' @importFrom ggplot2 ggplot
#' @importFrom rstatix ggscatterstats
#' @importFrom magrittr `%>%`
#' @importFrom reshape2 melt
#' @importFrom dplyr select
#' @importFrom cowplot ggdraw
#' @importFrom cowplot plot_grid
#'
#' @examples
#' getFinalClusters(SeuData)
getFinalClusters <- function(SeuData) {

  # Get the list of object names ending with _Clean or _Extra
  object_names <- ls(envir = globalenv())[grepl("_Clean$|_Extra$", ls(envir = globalenv()))]

  # Loop through the object names
  for (name in object_names) {
    # Extract the cluster name by removing the _Clean or _Extra suffix
    cluster_name <- gsub("(_Clean$)|(_Extra$)", "", name)

    # Get the cell barcodes from the data frame
    cell_barcodes <- row.names(get(name, envir = globalenv()))

    # Update the Idents of the clusters in the Seurat object
    SeuratObject::Idents(object = SeuData, cells = cell_barcodes) <- as.character(cluster_name)
  }

  # Update the celltype column in the Seurat object
  SeuData$celltype <- SeuratObject::Idents(SeuData)

  plot <- Seurat::DimPlot(SeuData, reduction = "umap", label = TRUE)
  print(plot)

  # Return the updated Seurat object
  return(SeuData)

  }
