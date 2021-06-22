#' Update LOC gene names based on Ensembl file
#'
#' @param matrix Matrix to fix
#' @param path Path of Ensembl MF5 file
#'
#' @return Matrix with new rownames
#' @export
UpdateGenes <- function(matrix, path = "data/reference/gene_info_MF5.txt") {

  ensembl_MF5 <- readr::read_tsv(path)
  ensembl_MF5 <- ensembl_MF5 %>% janitor::clean_names()
  genes <- rownames(matrix)
  LOCgenes <- data.frame(LOCID = genes[grep("^LOC", genes)])
  LOCgenes$ID <- gsub("LOC", "", LOCgenes$LOCID)
  LOCgenes <- merge(LOCgenes, ensembl_MF5, by.x = "ID", by.y = "ncbi_gene_formerly_entrezgene_id")
  LOCgenes_switch <- LOCgenes %>% filter(!is.na(gene_name)) %>% select(LOCID, gene_name) %>% distinct()
  matches <- match(LOCgenes_switch$LOCID, genes)
  rownames(matrix)[matches] <- as.character(LOCgenes_switch$gene_name)
  return(matrix)

}

#' Read counts from .txt file
#'
#' @param filename filename
#' @param gz TRUE/FALSE indicating file zip
#'
#' @return DGE matrix
#' @export
ReadTxtCounts <- function(filename, gz = TRUE) {
  if (gz == TRUE) dge <- readr::read_delim(gzfile(filename), delim = "\t", col_types = readr::cols(.default = readr::col_integer(), GENE = readr::col_character()))
  if (gz == FALSE) dge <- readr::read_delim(filename, delim = "\t", col_types = readr::cols(.default = readr::col_integer(), GENE = readr::col_character()))
  genes <- dge[, 1, drop = TRUE]
  barcodes <- colnames(dge)[2:ncol(dge)]
  dge <- dge[, -1]
  dge <- as.matrix(dge)
  mode(dge) <- "integer"
  dge <- Matrix::Matrix(dge, sparse = TRUE)
  rownames(dge) <- genes
  colnames(dge) <- barcodes
  return(dge)
}

RunSoupX <- function(
  seurat.obj,
  droplets,
  sample,
  cluster.var,
  plot = TRUE,
  save = FALSE
) {

  # check column names
  assertthat::assert_that(all(table(colnames(seurat.obj) %in% colnames(droplets))))
  seurat.obj <- AddMetaData(seurat.obj, metadata = as.data.frame(seurat.obj@reductions$pca_umap@cell.embeddings))
  toc <- seurat.obj@assays$RNA@counts
  dim(toc)
  colnames(toc)

  # pull droplet barcodes
  lib_sizes <- Matrix::colSums(droplets)
  colnames(toc) %in% colnames(droplets)
  keep_barcodes <- names(lib_sizes)[lib_sizes <= 100 | names(lib_sizes) %in% colnames(toc)]
  tod <- droplets[, keep_barcodes]
  dim(tod)

  # intersect genes
  gene_isx <- intersect(rownames(toc), rownames(tod))
  tod <- as.matrix(tod[gene_isx, ])
  toc <- as.matrix(toc[gene_isx, ])
  assertthat::assert_that(dim(tod)[1] == dim(toc)[1] & dim(tod)[2] > dim(toc)[2])

  # generate soup channel and set clusters
  sc <- SoupX::SoupChannel(tod, toc, metaData = seurat.obj@meta.data)
  sc <- SoupX::setClusters(sc, setNames(as.character(seurat.obj@meta.data[[cluster.var]]), colnames(seurat.obj)))
  sc <- SoupX::setDR(sc, seurat.obj@meta.data[colnames(sc$toc), c("pUMAP_1", "pUMAP_2")])

  # estimate contamination and generate adjusted counts
  sc <- SoupX::autoEstCont(sc)
  # adj_counts <- SoupX::adjustCounts(sc, roundToInt = TRUE)
  # cntSoggy = rowSums(sc$toc > 0)
  # cntStrained = rowSums(adj_counts > 0)
  # mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
  # tail(sort(rowSums(sc$toc > adj_counts)/rowSums(sc$toc > 0)), n = 20)
  # SoupX::plotChangeMap(sc, adj_counts, "CPA3")

  # generate figures and save objects
  if (plot) {
    plot_df <- sc$fit$dd
    plot_df <- plot_df %>% group_by(gene) %>% summarize(meanrho = mean(rhoEst), soupExp = mean(soupExp))
    nolist <- as.character(plot_df$gene[plot_df$meanrho >= quantile(plot_df$meanrho, 0.90)])
    topexpr <- head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 50)
    plot_df$label <- FALSE
    plot_df$label[plot_df$gene %in% nolist] <- TRUE
    a <- ggplot(plot_df, aes(x = soupExp, y = meanrho, fill = label)) +
      geom_point(shape = 21, color = "black", size = 3, alpha = 0.6) +
      scale_fill_manual(values = c("gray80", "#D55E00")) +
      guides(fill = FALSE) +
      ggrepel::geom_text_repel(data = plot_df[plot_df$label == TRUE, ], aes(label = gene),
                               size = 4, fontface = "bold", nudge_x = 5E-5, nudge_y = 1) +
      labs(y = "Ambient RNA % average per gene", x = "Soup expression", title = "Average gene contamination and soup expression") +
      GeneralTheme(base_size = 14)

    d <- ggplot(sc$fit$dd, aes(x = as.factor(rhoIdx), y = rhoEst)) +
      geom_boxplot() +
      labs(x = "Cluster", y = "Estimated contamination %", title = "Ambient RNA % across cluster") +
      GeneralTheme(base_size = 14)

    b <- ggplot(plot_df, aes(x = meanrho, y = ..scaled..)) +
      geom_density() +
      scale_y_continuous(breaks = c(0, 1)) +
      labs(x = "", y = "Density") +
      coord_flip() +
      GeneralTheme(base_size = 14)

    c <- DimPlot(obj, pt.size = 1, reduction = "umap", group.by = "base_batch_leiden",
                 label = TRUE, repel = TRUE, label.size = 6) +
      GeneralTheme(14) +
      labs(x = "UMAP1", y = "UMAP2", title = "UMAP embedding colored by cluster") +
      NoLegend()
    c$layers[[2]]$aes_params$fontface <- "bold"

    contam_grid <- a + b + c + d + plot_layout(nrow = 2, widths = c(1, 1))
    contam_grid
    ggsave(filename = glue("plots/{sample}_contamination_plot.png"), plot = contam_grid,
           width = 10, height = 10, units = "in", dpi = 300, scale = 1)
  }

  if (save) {
    #saveRDS(sc, file = glue("data/soupx/{sample}_soupx_object.rds"))
    #saveRDS(adj_counts, file = glue("data/soupx/{sample}_adj_counts.rds"))
  }

  # return objects
  fit <- sc
  fit$sample <- sample
  return(fit)
}
