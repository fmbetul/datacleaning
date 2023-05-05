global_vars <- ls(envir = globalenv())
extra_vars <- global_vars[grepl("_Extra$", global_vars)]
all_vars <- c("NewLeftovers", extra_vars)
utils::globalVariables(all_vars)
rm(global_vars)
rm(extra_vars)
rm(all_vars)

#' Assign cells in the _Extra cluster based on their UMAP coordinates
#'
#' This function assigns cells in the _Extra cluster based on their UMAP coordinates.
#'
#' @param SeuData Seurat object with UMAP coordinates
#' @param ClusterName Name of the cluster to assign the _Extra cells to
#' @param min.x Minimum x value for the _Extra cluster boundary
#' @param max.x Maximum x value for the _Extra cluster boundary
#' @param min.y Minimum y value for the _Extra cluster boundary
#' @param max.y Maximum y value for the _Extra cluster boundary
#' @return A Seurat object with the _Extra cluster assigned
#'
#' @importFrom Seurat DimPlot
#' @importFrom ggplot2 ggplot
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' assignExtraCluster(SeuratObject, ClusterName = "Cluster1", -10, 10, -10, 10)
assignExtraCluster <- function(SeuData, ClusterName, min.x = -100, max.x = 100, min.y = -100, max.y = 100) {

  # Get the UMAP coordinates of each cell
  UMAP_Data <- as.data.frame(SeuData@reductions$umap@cell.embeddings)

  # Check Leftovers
  Leftovers <- get("Leftovers", envir = globalenv())
  print(paste("Number of cells in current Leftovers:", dim(Leftovers)[1], sep = " "))
  print(CheckUMAP3(sample = SeuData, CheckInput = Leftovers))

  # Define the borders of the _Extra cluster
  ExtraCluster <- subset(UMAP_Data, row.names(UMAP_Data) %in% row.names(Leftovers)
                         & UMAP_Data$UMAP_1 > min.x
                         & UMAP_Data$UMAP_1 < max.x
                         & UMAP_Data$UMAP_2 > min.y
                         & UMAP_Data$UMAP_2 < max.y)

  # Visualize the _Extra cluster
  print(paste("Number of cells assigned to", ClusterName, "_Extra cluster:", dim(ExtraCluster)[1], sep = " "))
  print(CheckUMAP1(sample = SeuData, CheckInput = ExtraCluster))

  # Save the _Extra cluster
  set.name <- paste0(ClusterName, "_Extra")
  assign(set.name, ExtraCluster, envir = globalenv())

  # Visualize and save NewLeftovers in the global environment
  NewLeftovers <- subset(Leftovers, !(row.names(Leftovers) %in% row.names(ExtraCluster)))
  assign("NewLeftovers", NewLeftovers, envir = globalenv())
  print(paste("Number of cells in NewLeftovers must be:", dim(Leftovers)[1] - dim(ExtraCluster)[1], sep = " "))
  print(paste("Number of cells in NewLeftovers is:", dim(NewLeftovers)[1], sep = " "))
  print(CheckUMAP4(sample = SeuData, CheckInput = NewLeftovers))

  # look at all clusters overall
  print(Seurat::DimPlot(SeuData, reduction = "umap", label  =TRUE))

  }
