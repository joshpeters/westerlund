
<!-- README.md is generated from README.Rmd. Please edit that file -->

# westerlund

<!-- badges: start -->
<!-- badges: end -->

`westerlund` is a toolkit for single-cell RNA-sequencing analysis built
on top of Seurat providing functions to quickly and comprehensively
assess clustering results. Inspired by many publications and prior
packages (highlighted below), `westerlund` enables data sketching,
cluster bootstrapping, automatic and assisted resolution scanning, and
metric visualizations.

## Why `westerlund`?

The functions within `westerlund` provide additional critical
functionality for common single-cell analyses that have been implemented
in a range of outstanding publications, including:  
1. Supporting functions based on the recommendations by the authors of
[pipeComp](https://doi.org/10.1186/s13059-020-02136-7),  
1. Leiden clustering implemented in the `leidenbase` algorithm by the
[Trapnell Lab](https://github.com/cole-trapnell-lab/leidenbase),  
1. Interpretable resolution and k-nearest neighbors scanning for
determining the optimal cluster granularity utilizing the `future`
framework,  
1. Heuristic metrics to determine appropriate PCs to utilize for further
dimensionality reduction and clustering,  
1. Cluster robustness measures including SVM classification and adjusted
Rand index random seed and subset sampling,  
1. Functions to quickly assess cluster significance based on DE genes
and merge clusters based on size and/or user flagging,  
1. Gene list scoring using scaled expression instead of normalized
expression as utilized by [Smillie et
al.](https://doi.org/10.1016/j.cell.2019.06.029),  
1. Bayesian classifiers as implemented by [Zeymour et
al.](https://doi.org/10.1038/s41590-018-0051-0),  
1. Gene specificity index calculations as proposed by [Tosches et
al.](https://doi.org/10.1126/science.aar4237). 1. Max-to-second
expression ratios as utilized in [Zilions et
al.](https://doi.org/10.1016/j.immuni.2019.03.009)

These functions are plug and play with `Seurat` objects and can enable
easy expansion of available scRNA-seq analytical tools.

## Installation

    # install.packages("devtools")
    devtools::install_github("joshpeters/westerlund")

## Examples

``` r
library(tidyverse)
library(Seurat)
library(SeuratData)
library(westerlund)

#SeuratData::InstallData("pbmc3k")
#SeuratData::InstalledData()
# object <- SeuratData::LoadData("pbmc3k")
# object$seurat_annotations <- as.character(object$seurat_annotations)
# object$seurat_annotations[is.na(object$seurat_annotations)] <- "Unassigned"
# Idents(object) <- "seurat_annotations"
# object <- ReduceDims(seurat.obj = object, n.var.features = 3000, n.pcs = 50, remove.genes = NULL)
# object <- Seurat::RunUMAP(object = object, reduction = "pca", reduction.name = "pca_umap",
#                   dims = 1:object@misc$pca_metrics$tp_pc, seed.use = 1, verbose = FALSE)
# plot <- DimPlot(object, reduction = "pca_umap", group.by = "seurat_annotations", pt.size = 1, shuffle = TRUE, label = FALSE, label.size = 6, repel = TRUE)
# plot + theme_classic(base_size = 18) +
#   #NoLegend() +
#   RemoveAxes() + RemoveBackgrounds(outline = FALSE) + SpaceAxesTitles() +
#   colorspace::scale_color_discrete_qualitative("Dark 3") +
#   labs(x = "UMAP1", y = "UMAP2", title = "")
# markers <- DefineMarkers(object, "seurat_annotations")
# Bcell_markers <- markers %>% filter(group == "B") %>% top_n(100, auc) %>% pull(feature)
# object <- AddModuleScoreScaled(object, list(Bcell_markers), name = "SMS")
# object@meta.data[, grep("^MS", colnames(object[[]]))] <- NULL
# object <- AddModuleScore(object, list(Bcell_markers), name = "MS", nbin = 20, ctrl = 100)
```
