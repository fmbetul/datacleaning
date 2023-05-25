#' CheckUMAP function
#'
#' This function adds metadata to the input Seurat object, assigns the metadata to sample's Idents,
#' and then generates a UMAP plot of the Seurat object with highlighted cells of selected cluster.
#'
#' @param CheckInput A data.frame or object to be checked. If it is a data.frame,
#'                   the function generates metadata based on the row names.
#'                   If it is not a data.frame, the function generates metadata
#'                   based on the column names.
#' @param sample A Seurat object with UMAP coordinates. The default is `SeuData`.
#' @param color The color for the UMAP plot. The default is `#cb181d`.
#' @param label.text The text to use in generating the metadata. The default is `"RoI"`.
#'
#' @return A UMAP plot for the specified subset of cells
#' @importFrom Seurat DimPlot
#' @importFrom SeuratObject Idents AddMetaData
#'


CheckUMAP = function(CheckInput, sample, color = "#cb181d", label.text = "RoI"){

  if(inherits(CheckInput, "data.frame")){

    CheckMeta = as.data.frame(c(rep(label.text, length(row.names(CheckInput)))))
    colnames(CheckMeta) = "Pop"
    row.names(CheckMeta) = row.names(CheckInput)

  }else{

    CheckMeta = as.data.frame(c(rep(label.text, length(colnames(CheckInput)))))
    colnames(CheckMeta) = "Pop"
    row.names(CheckMeta) = colnames(CheckInput)

  }

  sample = SeuratObject::AddMetaData(sample, CheckMeta, "CheckMeta")
  SeuratObject::Idents(sample) = "CheckMeta"
  Seurat::DimPlot(sample, reduction="umap", cols = color, na.value = "grey85")

}
