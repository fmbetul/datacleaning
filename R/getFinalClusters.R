#' getFinalClusters function
#'
#' This function updates the identities of the cells in the input Seurat object
#' according to the "_Clean" and "_Extra" clusters defined in the global "datacleaning_list",
#' and plots the final UMAP visualization. It also updates the "celltype" column in
#' the Seurat object's metadata and removes the "datacleaning_list" from the global environment.
#'
#' @param SeuData A Seurat object that needs to have its cell identities updated.
#'
#' @return SeuData The updated Seurat object, with the new cell identities reflected in the metadata.
#' It also prints a UMAP plot of the updated Seurat object.
#' @importFrom Seurat DimPlot
#' @importFrom SeuratObject Idents
#' @export
#'
#' @examples
#' # Assuming 'SeuratObject' is your Seurat object:
#' SeuratObject <- getFinalClusters(SeuratObject)


getFinalClusters <- function(SeuData) {

  # Check if "datacleaning_list" exists in the global environment
  if (!exists("datacleaning_list", envir = globalenv())) {
    stop("Error: No datacleaning_list found in global environment.")
  }

  # Get "datacleaning_list" from the global environment
  datacleaning_list <- get("datacleaning_list", envir = globalenv())

  # Get the object names ending with _Clean or _Extra from "datacleaning_list"
  object_names <- grep("_Clean$|_Extra$", names(datacleaning_list), value = TRUE)

  # Loop through the object names
  for (name in object_names) {
    # Extract the cluster name by removing the _Clean or _Extra suffix
    cluster_name <- gsub("(_Clean$)|(_Extra$)", "", name)

    # Get the cell barcodes from the data frame
    cell_barcodes <- row.names(datacleaning_list[[name]])

    # Update the Idents of the clusters in the Seurat object
    SeuratObject::Idents(object = SeuData, cells = cell_barcodes) <- as.character(cluster_name)
  }

  # Update the celltype column in the Seurat object metadata
  SeuData$celltype <- SeuratObject::Idents(SeuData)

  plot <- Seurat::DimPlot(SeuData, reduction = "umap", label = TRUE)
  print(plot)

  # remove datacleaning_list from the global environment
  rm(list = "datacleaning_list", envir = globalenv())

  # Return the updated Seurat object
  return(SeuData)

}
