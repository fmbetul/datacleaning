utils::globalVariables(c("Leftovers"))

#' Validate and visualize the updated Leftovers after assigning Extra cluster
#'
#' @param SeuData Seurat object containing the data to be analyzed
#'
#' @return A plot of the UMAP coordinates of the updated Leftovers
#'
#' @import Seurat
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_size_continuous
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_shape_manual
#'
#' @export
#'
#' @examples
#' validateExtraCluster(SeuratObject)
validateExtraCluster <- function(SeuData){

  # Update Leftovers
  assign("Leftovers", NewLeftovers, envir = globalenv())

  print(paste("Number of cells in updated Leftovers must be:", dim(NewLeftovers)[1], sep = " "))
  print(paste("Number of cells in updated Leftovers:", dim(Leftovers)[1], sep = " "))
  print(CheckUMAP3(sample = SeuData, CheckInput = Leftovers))

  }
