#' assignNewCluster function
#'
#' This function selects cells from the leftovers (cells not assigned to any cluster yet) based on
#' provided x and y coordinates, assigns them to a new cluster, and updates the leftovers.
#' In order for the new cells to be successfully assigned and subtracted from the unassigned cells,
#' the validate parameter must be set to TRUE.
#'
#' @param SeuData A Seurat object that contains the data of the clusters.
#' @param ClusterName Name of the new cluster to be generated.
#' @param min.x Minimum x coordinate of the new cluster boundaries. The default is `-100`.
#' @param max.x Maximum x coordinate of the new cluster boundaries. The default is `100`.
#' @param min.y Minimum y coordinate of the new cluster boundaries. The default is `-100`.
#' @param max.y Maximum y coordinate of the new cluster boundaries. The default is `100`.
#' @param validate Boolean value indicating whether to validate the assignment. This must be set to `TRUE`
#'                 for the new cells to be successfully assigned and subtracted from the unassigned cells.
#'                 The default is `FALSE`.
#'
#' @return Invisibly returns a Seurat UMAP plot of all clusters in the data. This function also updates the
#'         global list "datacleaning_list" by adding a new cluster and updating the leftovers.
#' @importFrom Seurat DimPlot
#' @export
#'


assignNewCluster <- function(SeuData, ClusterName, min.x = -100, max.x = 100, min.y = -100, max.y = 100, validate = FALSE) {


  # Check if "datacleaning_list" exists in the global environment
  if (!exists("datacleaning_list", envir = globalenv())) {
    stop("Error: No datacleaning_list found in global environment.")
  }

  # Get "datacleaning_list" from the global environment
  datacleaning_list <- get("datacleaning_list", envir = globalenv())


  # Check if "Leftovers" exists in "datacleaning_list"
  if (!"Leftovers" %in% names(datacleaning_list)) {

    # Get all _Leftover objects from datacleaning_list
    leftover_objects <- grep("_Leftover", names(datacleaning_list), value = TRUE)

    # If there are no _Leftover objects in datacleaning_list, stop the function
    if (length(leftover_objects) == 0) {
      stop("Error: No _Leftover objects found in datacleaning_list. Crop at least one cluster using cropCluster function.")
    }

    # Merge all _Leftover objects to generate a single Leftovers object
    Leftovers <- do.call(rbind, datacleaning_list[leftover_objects])

    # Remove the prefix from the row names
    rownames(Leftovers) <- sub(".*\\.", "", rownames(Leftovers))

    # Save Leftovers object to datacleaning_list in the global environment
    datacleaning_list[["Leftovers"]] <- Leftovers
    assign("datacleaning_list", datacleaning_list, envir = globalenv())
  }

  # Get "Leftovers" from "datacleaning_list"
  Leftovers <- datacleaning_list[["Leftovers"]]

  # Subset "Leftovers" based on provided x and y boundaries
  ExtraCluster <- subset(Leftovers,
                         Leftovers$UMAP_1 > min.x & Leftovers$UMAP_1 < max.x &
                           Leftovers$UMAP_2 > min.y & Leftovers$UMAP_2 < max.y)

  # Save _Extra cluster to "datacleaning_list" in the global environment
  set.extra.name <- paste0(ClusterName, "_Extra")
  datacleaning_list[[set.extra.name]] <- ExtraCluster
  assign("datacleaning_list", datacleaning_list, envir = globalenv())

  # Calculate NewLeftovers and save to "datacleaning_list" in the global environment
  NewLeftovers <- subset(Leftovers, !(rownames(Leftovers) %in% rownames(ExtraCluster)))
  datacleaning_list[["NewLeftovers"]] <- NewLeftovers
  assign("datacleaning_list", datacleaning_list, envir = globalenv())


  # output
  print(paste("Number of cells in current Leftovers:", dim(Leftovers)[1], sep = " "))
  print(CheckUMAP(sample = SeuData, CheckInput = Leftovers, color = "purple", label.text = "unassigned"))
  # ExtraCluster
  print(paste("Number of cells assigned to new cluster:", dim(ExtraCluster)[1], sep = " "))
  print(CheckUMAP(sample = SeuData, CheckInput = ExtraCluster, color = "#008000", label.text = ClusterName))
  # NewLeftovers
  print(paste("Number of cells in NewLeftovers is:", dim(NewLeftovers)[1], sep = " "))
  print(CheckUMAP(sample = SeuData, CheckInput = NewLeftovers, color = "#cb181d", label.text = "new_unassigned"))
  # DimPlot
  print(Seurat::DimPlot(SeuData, reduction = "umap", label  =TRUE))

  if(validate == TRUE){

    # Get "datacleaning_list" from the global environment
    datacleaning_list <- get("datacleaning_list", envir = globalenv())

    # Get "NewLeftovers" from "datacleaning_list"
    NewLeftovers <- datacleaning_list[["NewLeftovers"]]

    # Save NewLeftovers as "Leftovers" to "datacleaning_list" in the global environment
    datacleaning_list[["Leftovers"]] <- NewLeftovers
    assign("datacleaning_list", datacleaning_list, envir = globalenv())

    # output
    print("Cluster validated ! Leftovers updated ...")
    print(paste("Number of cells in updated Leftovers:", dim(datacleaning_list[["Leftovers"]])[1], sep = " "))

  }else{

    # remove _Extra cluster from the datacleaning_list
    datacleaning_list <- get("datacleaning_list", envir = globalenv())  # get the current "datacleaning_list"
    datacleaning_list[[set.extra.name]] <- NULL  # remove the item from the list
    assign("datacleaning_list", datacleaning_list, envir = globalenv())  # update "datacleaning_list" in the global environment

  }

}
