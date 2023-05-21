#' cropCluster function
#'
#' This function selects a cluster from the Seurat object based on the provided cluster name or number,
#' adjusts its boundaries based on provided min and max coordinates in both axes, and then stores
#' the cleaned cluster and leftovers in datacleaning_list in the global environment.
#'
#' @param SeuData A Seurat object that contains the data of the clusters.
#' @param ClusterName Name of the cluster to be selected.
#' @param ClusterNo Number of the cluster to be selected. The default is `0`.
#' @param min.x Minimum x coordinate of the new cluster boundaries. The default is `-100`.
#' @param max.x Maximum x coordinate of the new cluster boundaries. The default is `100`.
#' @param min.y Minimum y coordinate of the new cluster boundaries. The default is `-100`.
#' @param max.y Maximum y coordinate of the new cluster boundaries. The default is `100`.
#'
#' @return Invisibly returns a Seurat UMAP plot of all clusters in the data.
#'        This function also generates two global variables, one for the cleaned cluster
#'        and one for the leftover cells, and adds them to the global list "datacleaning_list".
#' @importFrom Seurat DimPlot
#' @importFrom SeuratObject Idents WhichCells
#' @export
#'


cropCluster = function(SeuData, ClusterName, ClusterNo = 0,
                       min.x = -100, max.x = 100,
                       min.y = -100, max.y = 100){



  # if "datacleaning_list" does not exist in the global environment, generate it
  if (!exists("datacleaning_list", envir = globalenv())) {
    assign("datacleaning_list", list(), envir = globalenv())
  }

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
  check.cluster <- CheckUMAP(Cluster_Coords, sample = SeuData, color = "blue", label.text = ClusterName)
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
  print(paste("Number of cells in Current_Leftovers will be:", length(cluster.barcodes) - dim(Cluster_Clean)[1], sep = " "))

  # look at the new cluster on UMAP
  check.cropped.cluster <- CheckUMAP(Cluster_Clean, sample = SeuData, color = "#008000", label.text = paste0("new_", ClusterName))
  print(check.cropped.cluster)

  # Save the _Clean cluster to the "datacleaning_list"
  set.clean.name <- paste0(ClusterName, "_Clean")     # eg."POMC_Clean"

  # get the current "datacleaning_list"
  datacleaning_list <- get("datacleaning_list", envir = globalenv())

  # add the new object to the list
  datacleaning_list[[set.clean.name]] <- Cluster_Clean

  # update "datacleaning_list" in the global environment
  assign("datacleaning_list", datacleaning_list, envir = globalenv())


  ###

  # 3. Current_Leftovers (double-check)

  Current_Leftovers <- subset(Cluster_Coords, !(row.names(Cluster_Coords) %in% row.names(Cluster_Clean)))
  print(paste("Current number of cells in current Leftovers is", dim(Current_Leftovers)[1], sep = " "))

  # look at the Leftovers on UMAP
  check.leftovers <- CheckUMAP(Current_Leftovers, sample = SeuData, color = "#cb181d", label.text = "removed_cells")
  print(check.leftovers)

  # Save the _Leftover cluster to the "datacleaning_list"
  set.leftover.name <- paste0(ClusterName, "_Leftover")     # eg."POMC_Leftover"

  # get the current "datacleaning_list"
  datacleaning_list <- get("datacleaning_list", envir = globalenv())

  # add the new object to the list
  datacleaning_list[[set.leftover.name]] <- Current_Leftovers

  # update "datacleaning_list" in the global environment
  assign("datacleaning_list", datacleaning_list, envir = globalenv())


  # look at all clusters overall
  print(Seurat::DimPlot(SeuData, reduction = "umap", label  =TRUE))

}
