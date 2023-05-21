# datacleaning

The `datacleaning` package offers a suite of functions designed to enhance the process of cleaning up and visualizing clusters within a `Seurat` object, commonly used for single-cell RNA sequencing (scRNA-seq) data analysis. It includes functions for cropping clusters, reassigning cells to new clusters based on UMAP coordinates, visualizing the refined clusters, updating cell identities based on specified clusters and visualizing the resulting data using UMAP plots.


## Installation

To install the package from GitHub, utilize the `devtools` package:

``` r
devtools::install_github("fmbetul/datacleaning")
```

## Usage

To use the package, load it into your R session:

``` r
library(datacleaning)
```

Below is an example showcasing the use of various functions within this package:

``` r
# Assuming `SeuratObject` is your Seurat object
# Step 1: Apply CropCluster() to each cluster needing cropping, either by providing ClusterNo or ClusterName
cropCluster(SeuratObject, ClusterName = "Cluster1", min.x = -10, max.y = 10)
cropCluster(SeuratObject, ClusterName = "Cluster2", min.x = -5, min.y = 5)
```

```r
# Step 2: Assign leftover cells to a new cluster, set validate = TRUE after you feel satisfied with the new cluster.
assignNewCluster(SeuratObject, ClusterName = "Cluster3", min.x = -10, max.y = 10, validate = TRUE)
assignNewCluster(SeuratObject, ClusterName = "Cluster4", min.x = -5, min.y = -5, validate = TRUE)
```

```r
# Step 3: Obtain the final cleaned-up clusters
SeuratObject <- getFinalClusters(SeuratObject)
```

For more information on available functions and their usage, please refer to the package documentation.


## Dependencies

The `datacleaning` package has the following dependencies:

- R (>= 3.6.0)
- Seurat (>= 4.0.0)


## License

This package is licensed under the GNU General Public License version 3 (GPL-3). See the [LICENSE](LICENSE) file for more details.


## Issues and Contributions

If you encounter any issues with the `datacleaning` package or would like to contribute to its development, please [open an issue](https://github.com/fmbetul/datacleaning/issues) on GitHub.


## Contact

For any further inquiries or questions, please feel free to contact the package maintainer, F. M. Betul Erol, at fmab_erol@hotmail.com.
