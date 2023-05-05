#' Check UMAP plot for a specific subset of cells
#'
#' This function checks the UMAP plot for a specific subset of cells.
#'
#' @param CheckInput A data frame or matrix with UMAP coordinates to check
#' @param sample A Seurat object with UMAP coordinates
#' @return A UMAP plot for the specified subset of cells
#'
#' @importFrom Seurat DimPlot
#' @importFrom SeuratObject Idents
#' @importFrom ggplot2 ggplot
#' @importFrom magrittr %>%
#'
CheckUMAP1 = function(CheckInput, sample = SeuData){

  if(inherits(CheckInput, "data.frame")){

    CheckMeta = as.data.frame(c(rep("current_cluster", length(row.names(CheckInput)))))
    colnames(CheckMeta) = "Pop"
    row.names(CheckMeta) = row.names(CheckInput)

  }else{

    CheckMeta = as.data.frame(c(rep("current_cluster", length(colnames(CheckInput)))))
    colnames(CheckMeta) = "Pop"
    row.names(CheckMeta) = colnames(CheckInput)
    }

  sample = SeuratObject::AddMetaData(sample, CheckMeta, "CheckMeta")
  SeuratObject::Idents(sample) = "CheckMeta"
  Seurat::DimPlot(sample, reduction="umap", cols = c("#cb181d"), na.value = "grey85")

  }
