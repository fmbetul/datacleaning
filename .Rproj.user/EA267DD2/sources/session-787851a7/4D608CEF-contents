global_vars <- ls(envir = globalenv())
clean_vars <- global_vars[grepl("_Clean$", global_vars)]
all_vars <- c("Leftovers", clean_vars)
utils::globalVariables(all_vars)
rm(global_vars)
rm(clean_vars)
rm(all_vars)

#' Crop a cluster based on its UMAP coordinates
#'
#' This function crops a cluster based on its UMAP coordinates and saves the cleaned cluster as a new cluster with a specified name.
#'
#' @param SeuData A Seurat object containing the data to be cropped.
#' @param ClusterName A character string specifying the name of the cluster to be cropped.
#' @param ClusterNo An integer specifying the number of the cluster to be cropped. If specified, the function will ignore the ClusterName parameter and use the ClusterNo instead. Default is 0.
#' @param min.x A numeric value specifying the minimum x-coordinate of the cropped cluster on the UMAP plot. Default is -100.
#' @param max.x A numeric value specifying the maximum x-coordinate of the cropped cluster on the UMAP plot. Default is 100.
#' @param min.y A numeric value specifying the minimum y-coordinate of the cropped cluster on the UMAP plot. Default is -100.
#' @param max.y A numeric value specifying the maximum y-coordinate of the cropped cluster on the UMAP plot. Default is 100.
#'
#' @return Nothing. Saves the cleaned cluster as a new cluster with a specified name.
#'
#' @import Seurat
#' @importFrom SeuratObject Idents
#' @importFrom ggplot2 ggplot
#' @importFrom magrittr %>%
#'
#' @examples
#' CropCluster(SeuData, ClusterName = "Cluster1", min.x = -50, max.x = 50, min.y = -50, max.y = 50)
#'
#' @export
CropCluster = function(SeuData, ClusterName, ClusterNo = 0,
                       min.x = -100, max.x = 100,
                       min.y = -100, max.y = 100){


  # get cluster names:
  Clusters <- levels(SeuratObject::Idents(SeuData))

  # get the UMAP coordinates of each cell
  UMAP_Data <- as.data.frame(SeuData@reductions$umap@cell.embeddings)




  # 1. Select a new cluster & check its properties

  if(ClusterNo > 0){
    ClusterName <- Clusters[ClusterNo]
    if (!is.character(ClusterName)) {
      stop("Error: This ClusterNo is not available.")
    }
  }

  if(!(ClusterName %in% Clusters)){
    stop("Error: This ClusterName is not valid.")
  }

  cluster.barcodes <- SeuratObject::WhichCells(SeuData, idents = ClusterName)

  # Check-point
  print(paste("Current cluster:", ClusterName, sep = " "))
  print(paste("Total number of cells in this cluster:", length(cluster.barcodes), sep = " "))

  # get UMAP coordinates of the cells in the selected cluster
  Cluster_Coords <- subset(UMAP_Data, row.names(UMAP_Data) %in% cluster.barcodes)

  # look at the selected cluster on UMAP
  check.cluster <- CheckUMAP1(Cluster_Coords, sample = SeuData)
  print(check.cluster)

  ###

  # set the cluster boundaries
  Cluster_Clean = subset(UMAP_Data,
                         row.names(UMAP_Data) %in% cluster.barcodes
                         & UMAP_Data $UMAP_1 > min.x
                         & UMAP_Data $UMAP_1 < max.x
                         & UMAP_Data $UMAP_2 > min.y
                         & UMAP_Data $UMAP_2 < max.y)


  print(paste("Cells left in", ClusterName ,"cluster:", dim(Cluster_Clean)[1], sep = " "))
  print(paste("Current number of cells in current Leftovers must be", length(cluster.barcodes) - dim(Cluster_Clean)[1], sep = " "))

  # look at the new cluster on UMAP
  check.cropped.cluster <- CheckUMAP2(Cluster_Clean, sample = SeuData)
  print(check.cropped.cluster)

  # Save the _Clean cluster
  set.name <- paste0(ClusterName, "_Clean")     # eg."POMC_Clean"
  assign(set.name, Cluster_Clean, envir = globalenv())


  ###

  # 3. Current_Leftovers (double-check)

  Current_Leftovers <- subset(Cluster_Coords, !(row.names(Cluster_Coords) %in% row.names(Cluster_Clean)))
  print(paste("Current number of cells in current Leftovers is", dim(Current_Leftovers)[1], sep = " "))

  # look at the Leftovers on UMAP
  check.leftovers <- CheckUMAP3(Current_Leftovers, sample = SeuData)
  print(check.leftovers)

  # look at all clusters overall
  print(Seurat::DimPlot(SeuData, reduction = "umap", label  =TRUE))

  }
