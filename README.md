# datacleaning: Streamline Single-cell RNA-seq Data Clustering with Seurat

The `datacleaning` package offers a suite of functions designed to enhance the process of cleaning up and visualizing clusters within a `Seurat` object, which is generated from single-cell RNA-seq data. By incorporating functionalities for cropping clusters, reassigning cells to new clusters based on UMAP coordinates, and visualizing the refined clusters, this package streamlines the data analysis process.

## Installation

To install the package from GitHub, utilize the `devtools` package:

``` r
devtools::install_github("fmbetul/datacleaning")
```

## Usage

Below is an example showcasing the use of various functions within this package:

``` r
# Load required packages and data
library(Seurat)
library(datacleaning)

# Step 1: Crop clusters (apply CropCluster() to each cluster needing cropping, either by providing ClusterNo or ClusterName)
CropCluster(SeuratObject, ClusterNo = 1, min.x = -1, max.x = 4, max.y = -2)
CropCluster(SeuratObject, ClusterName = "Cluster1", min.x = -1, max.x = 4, max.y = -2)

# Step 2: Visualize remaining cells in cropped clusters
getLeftovers(SeuratObject)

# Step 3: Assign leftover cells to a cluster (sequentially apply assignExtraCluster() and validateExtraCluster() functions)
assignExtraCluster(SeuratObject, "Cluster1", min.x = 1.5, max.y = 1.4)
validateExtraCluster(SeuratObject)
assignExtraCluster(SeuratObject, "Cluster2", min.x = -3.1, min.y = 1.4)
validateExtraCluster(SeuratObject)

# Step 4: Obtain the final cleaned-up clusters
SeuratObject <- getFinalClusters(SeuratObject)
```

For further details on available functions, including input and output specifications, please refer to the package documentation.

## License

The `datacleaning` package is licensed under the GPL-3 license.
