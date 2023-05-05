utils::globalVariables(c("Leftovers"))

#' Get Leftover cells not assigned to any clusters
#'
#' This function gets the UMAP coordinates of each cell and retrieves the cell barcodes in the Seurat object that are not assigned to any clusters. It then saves the Leftovers object in the global environment and returns the result.
#'
#' @param SeuData A Seurat object containing the UMAP coordinates of each cell.
#'
#' @return A data frame containing the UMAP coordinates of each leftover cell.
#'
#' @export
#'
#' @importFrom Seurat WhichCells DimPlot
#' @importFrom ggplot2 ggplot
#' @importFrom magrittr %>%
#'
#' @examples
#' getLeftovers(SeuratObject)
getLeftovers <- function(SeuData){

  # get the UMAP coordinates of each cell
  UMAP_Data <- as.data.frame(SeuData@reductions$umap@cell.embeddings)

  # get the vector containing names of the _Clean objects
  Clean_Cluster_Names <- ls(envir = globalenv())[grepl("_Clean$", ls(envir = globalenv()))]

  # get the cell barcodes in the _Clean objects
  assigned_barcodes <- character(0)
  for(n in Clean_Cluster_Names){
    assigned_barcodes <- c(assigned_barcodes, row.names(get(n, envir = globalenv())))
  }

  print(paste("Number of cells assigned:", length(assigned_barcodes), sep = " "))
  print(paste("Number of cells in Leftovers must be:", dim(UMAP_Data)[1] - length(assigned_barcodes), sep = " "))

  # get the cell barcodes of the Leftovers
  Leftovers <- subset(UMAP_Data, !(row.names(UMAP_Data) %in% assigned_barcodes))

  # Save the Leftovers object in the global environment
  assign("Leftovers", Leftovers, envir = globalenv())

  print(paste("Number of cells in Leftovers:", dim(Leftovers)[1], sep = " "))
  print(CheckUMAP3(sample = SeuData, CheckInput = Leftovers))

  }
