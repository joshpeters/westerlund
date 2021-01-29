# FUNCTIONS ---------------------------------------------------------------

#' Leiden graph-based clustering
#'
#' leidenbase wrapper function to parse and select results
#'
#' @param res.param Resolution
#' @param graph Graph
#' @param seed.param Seed
#' @param return.params Return cluster-level results
#' @param return.clusters Return cell-level results
#' @param num.iter Number of leidenbase iterations
#'
#' @return List containing clustering and/or parameter results
#' @export
#'
#' @examples
Cluster <- function(res.param, graph, seed.param, return.params = TRUE, return.clusters = FALSE, num.iter = 1) {

  cluster_res <- leidenbase::leiden_find_partition(graph,
                                                   partition_type = "CPMVertexPartition",
                                                   initial_membership = NULL,
                                                   edge_weights = NULL,
                                                   node_sizes = NULL,
                                                   seed = seed.param,
                                                   resolution_parameter = res.param,
                                                   num_iter = num.iter,
                                                   verbose = FALSE)

  param_res <- data.frame(
    res_param = res.param,
    quality = cluster_res[["quality"]],
    modularity = cluster_res[["modularity"]],
    significance = cluster_res[["significance"]],
    cluster_count = max(cluster_res[["membership"]]))

  if(return.params & return.clusters) {
    return(list(cluster_res, param_res))
  } else if (return.params == TRUE & return.clusters == FALSE) {
    return(param_res)
  } else if (return.params == FALSE & return.clusters == TRUE) {
    return(cluster_res)
  }
}

ScanResolutions <- function(g,
                            res.start = 1E-10,
                            res.end = 1,
                            num.res = 100,
                            clusters.lb = 10,
                            clusters.ub = 50,
                            seed.param = 1,
                            use.mod = FALSE,
                            mod.percentile = 50) {

  res_params <-  signif(exp(seq(log(res.start), log(res.end), length.out = num.res)), 3)
  usethis::ui_todo("Scanning resolutions ({res_params[1]}, {res_params[length(res_params)]})")
  params_res <- future.apply::future_lapply(X = res_params, FUN = Cluster, graph = g, seed.param = seed.param, return.params = TRUE)
  params_res <- do.call(rbind, params_res)

  if (use.mod) {
    max_mod <- max(params_res$modularity)
    mod_threshold <- quantile(params_res$modularity[params_res$modularity > 0], mod.percentile)
    lower_mod_res <- min(params_res$res_param[params_res$modularity >= mod_threshold])
    upper_mod_res <- max(params_res$res_param[params_res$modularity >= mod_threshold])
    params <- c(lower_mod_res, upper_mod_res)
  } else {
    lower_cluster_res <- min(params_res$res_param[params_res$cluster_count >= clusters.lb])
    upper_cluster_res <- max(params_res$res_param[params_res$cluster_count <= clusters.ub])
    params <- c(lower_cluster_res, upper_cluster_res)
  }
  usethis::ui_done("Found bounds ({params[1]}, {params[2]})")
  return(params)
}

#' Create igraph-compatible graph and save in Seurat object
#'
#' @param object Seurat object
#' @param reduction Reduced dimension slot to pull from
#' @param dims Dimensions to use
#' @param knn Number of k-nearest neighbors to use
#'
#' @return Seurat object with igraph graph stored in `object@misc$westerlund_graph`
#' @export
PrepareGraph <- function(object, reduction, dims, knn) {
  graph.name <- glue::glue("{reduction}_snn_{knn}")
  stopifnot(ncol(object@reductions[[reduction]]) > dims)
  object <- Seurat::FindNeighbors(object = object, k.param = knn, prune.SNN = 1/15, dims = 1:dims,
                                  reduction = reduction, graph.name = c("nn", graph.name), compute.SNN = TRUE, verbose = FALSE)
  adj_matrix <- Matrix(as.matrix(object@graphs[[graph.name]]), sparse = TRUE)
  g <- graph_from_adjacency_matrix(adjmatrix = adj_matrix, mode = "undirected", weighted = TRUE, add.colnames = TRUE)
  object@misc$westerlund_graph <- g
  return(object)
}

SampleResolutions <- function(graph, res.start, res.end, num.res, num.samples) {

  res_params <-  signif(exp(seq(log(res.start), log(res.end), length.out = num.res)), 3)
  usethis::ui_todo("Sampling resolutions ({res_params[1]}, {res_params[length(res_params)]}), {length(rep(res_params, each = num.samples))} total")
  results <- future.apply::future_lapply(rep(res_params, each = num.samples), Cluster, graph = graph, seed.param = NULL,
                                         return.params = TRUE, return.clusters = TRUE, num.iter = 1)
  names(results) <- paste0("R_", rep(res_params, each = num.samples), "_", seq(1:num.samples))
  return(results)
}

VisualizeARI <- function(results) {

  res <- unique(stringr::str_match(names(results), "R_0\\.([[:digit:]]{2,8})_")[, 2, drop = TRUE])
  clusters <- lapply(results, `[[`, 1)
  aris <- map_df(.x = res, .f = ~ {
    mems <- lapply(clusters[grepl(.x, names(clusters))], `[[`, 1)
    base <- mems[[1]]
    aris <- sapply(mems[2:length(mems)], function(x, base) mclust::adjustedRandIndex(base, x), base = base)
    aris <- qdapTools::list2df(aris, col1 = "ari", col2 = "iter")
    aris$res <- .x
    return(aris)
  })
  a <- ggplot(aris, aes(x = ari, y = res)) +
    geom_jitter(alpha = 0.5, shape = 21, color = "black", fill = "white", position = position_jitter(height = 0.1)) +
    geom_boxplot(color = "black", fill = "transparent", outlier.shape = NA) +
    scale_x_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, 0.1)) +
    labs(x = "Adjusted Rand Index", y = "Resolution Parameter") +
    theme_classic(18) +
    theme(rec = element_rect(fill = "transparent"),
          axis.text = element_text(color = "black"),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 18),
          axis.title.y = element_text(margin = margin(0, 18, 0, 0)),
          axis.title.x = element_text(margin = margin(18, 0, 0, 0)),
          panel.grid.major.y = element_line(color = "gray50", size = 0.5, linetype = "dotted"),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent")
    )

  mods <- lapply(results, `[[`, 2)
  mods <- data.table::rbindlist(mods)
  mods$res_param <- as.factor(mods$res_param)
  b <- ggplot(mods, aes(x = modularity, y = res_param, group = res_param)) +
    geom_jitter(alpha = 0.5, shape = 21, color = "black", fill = "white", position = position_jitter(height = 0.2)) +
    geom_boxplot(color = "black", fill = "transparent", outlier.shape = NA) +
    #scale_x_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, 0.1)) +
    labs(x = "Modularity") +
    theme_classic(18) +
    theme(rec = element_rect(fill = "transparent"),
          axis.text = element_text(color = "black"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 18),
          axis.title.x = element_text(margin = margin(18, 0, 0, 0)),
          panel.grid.major.y = element_line(color = "gray50", size = 0.5, linetype = "dotted"),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent")
    )

  c <- ggplot(mods %>% group_by(res_param) %>% summarize(cluster_count = median(cluster_count)),
              aes(x = cluster_count, y = res_param, group = res_param)) +
    geom_col(color = "black", fill = "transparent") +
    geom_text(aes(x = cluster_count + 20, y = res_param, label = cluster_count), size = 4, fontface = "bold") +
    scale_x_continuous(limits = c(0, max(mods$cluster_count) + 30)) +
    labs(x = "# Clusters") +
    theme_classic(18) +
    theme(rec = element_rect(fill = "transparent"),
          axis.text = element_text(color = "black"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 18),
          axis.title.x = element_text(margin = margin(18, 0, 0, 0)),
          #panel.grid.major.y = element_line(color = "gray50", size = 0.5, linetype = "dotted"),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent")
    )
  grid <- a + plot_spacer() + b + plot_spacer() + c + plot_layout(widths = c(1, -0.1, 1, -0.1, 0.75))
  return(grid)
  #ggsave(filename = "~/bunchacruncha_test.pdf", plot = grid, bg = "transparent", width = 8, height = 4)


}

TestSignificance <- function(object, res, membership, thresholds = c(log(1.1), 0.05, 10)) {

  object$test_membership <- membership
  markers <- presto::wilcoxauc(object, group_by = "test_membership", assay = "data", seurat_assay = "RNA")
  group_statistics <- markers %>%
    filter(logFC > thresholds[1] & padj <= thresholds[2] & pct_in >= thresholds[3]) %>%
    group_by(group, .drop = FALSE) %>%
    summarize(num_DEGs = n()) %>%
    mutate(resolution = res)
  return(group_statistics)

}

WSS <- function(d) {
  sum(scale(d, scale = FALSE)^2)
}

Wrap <- function(membership, x) {
  spl <- split(x, membership)
  wss <- sum(sapply(spl, WSS))
}


PlotUMAP <- function(object, column, column.order, reduction = "umap") {

  object@meta.data <- object@meta.data %>%
    mutate(group = as.factor(as.character(!!sym(column))))
  object@meta.data$group <- fct_relevel(object@meta.data$group, column.order)
  object@meta.data$group <- droplevels(object@meta.data$group)
  object$umap1 <- object@reductions[[reduction]]@cell.embeddings[, 1, drop = TRUE]
  object$umap2 <- object@reductions[[reduction]]@cell.embeddings[, 2, drop = TRUE]
  centroids <- object@meta.data %>% arrange(group) %>%
    group_by(group) %>%
    summarize(umap1 = median(umap1), umap2 = median(umap2))
  centroids <- centroids %>% filter(group %in% column.order)

  plot <- ggplot(object@meta.data %>% arrange(desc(group)), aes(x = umap1, y = umap2, color = group)) +
    #ggrastr::geom_point_rast(alpha = 0.8, size = 2, raster.dpi = 300, color = "black") +
    ggrastr::geom_point_rast(size = 1, alpha = 0.5, raster.dpi = 300) +
    geom_point(centroids, mapping = aes(x = umap1, y = umap2, fill = group), size = 3, shape = 21, color = "black") +
    ggrepel::geom_text_repel(centroids, mapping = aes(x = umap1, y = umap2, label = group),
                             color = "black", size = 6,
                             max.overlaps = Inf, min.segment.length = unit(0.5, "lines"), point.padding = unit(2, "lines")) +
    labs(x = "UMAP1", y = "UMAP2") +
    colorspace::scale_color_discrete_qualitative("Dark 3") +
    colorspace::scale_fill_discrete_qualitative("Dark 3") +
    GeneralTheme(18) +
    theme(legend.position = "none",
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())
  #plot$layers[[3]]$aes_params$fontface <- "bold"
  return(plot)
}

GeneralTheme <- function(base_size, ...) {
  theme_classic(base_size = base_size, ...) +
    ggeasy::easy_all_text_color("black") +
    theme (
      axis.title.x = element_text(face = "plain", margin = margin(16, 0, 0, 0)),
      axis.title.y = element_text(face = "plain", margin = margin(0, 16, 0, 0), angle = 90, hjust = 0.5, vjust = 0.5),
      axis.line = element_blank(),
      plot.title = element_text(size =  base_size, color = "black", face = "bold", margin = margin(0,0,4,0)),
      plot.subtitle = element_text(size = base_size - 2, color = "black", margin = margin(0,0,4,0)),
      panel.background = element_rect(fill = "transparent", color = "black", size = 1),
      plot.background = element_rect(fill = "transparent", color = "transparent", size = 0),
      panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
      plot.caption = element_text(hjust = 0, color = "gray40", margin = margin(12)),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      legend.background = element_rect(fill = "transparent", color = "transparent"),
      legend.box.background = element_rect(fill = "transparent", color = "transparent"),
      legend.position = "right",
      legend.justification = "top",
      legend.key.size = unit(1, "line"),
      validate = TRUE
    )
}

ScanWalktrapKNN <- function(object, reduction = "pca", dims = 20, k.params = c(5, 10, 15, 20, 30, 50)) {

  # prepare graph
  usethis::ui_todo("Calculating {length(k.params)} graphs")
  graphs <- pbapply::pblapply(k.params, function(x, dims.param, reduction.param, object.param) {
    return(PrepareGraph(object = object.param, reduction = reduction.param, knn = x, dims = dims.param)@misc$westerlund_graph)
  }, dims.param = dims, reduction.param = reduction, object.param = object)
  barcodes <- Cells(object)

  # perform walktrap clustering
  usethis::ui_todo("Performing Walktrap clustering")
  walktraps <- pbapply::pblapply(graphs, function(x) {
    return(igraph::cluster_walktrap(x, steps = 4))
  })
  names(walktraps) <- as.character(k.params)
  object@misc$walktrap_results <- walktraps

  mems <- map_dfc(walktraps, membership)
  colnames(mems) <- paste0("walktrap_knn", as.character(round(1/k.params, 3)))
  rownames(mems) <- barcodes
  object@meta.data[, grep("walktrap", colnames(object[[]]))] <- NULL
  object <- AddMetaData(object, mems)
  return(object)

}

VisualizeSWK <- function(object) {

  stopifnot(!(is.null(object@misc$walktrap_results)))
  walktraps <- object@misc$walktrap_results
  k.params <- names(walktraps)

  sizes <- map_dfr(.x = k.params, .f = ~ {
    df <- data.frame(
      num_clusters = length(unique(walktraps[[.x]]$membership)),
      min_size = min(table(walktraps[[.x]]$membership)),
      max_size = max(table(walktraps[[.x]]$membership)),
      modularity = modularity(walktraps[[.x]]),
      knn = .x)
    return(df)
  })

  a <- ggplot(sizes, aes(x = num_clusters, y = modularity, label = knn)) +
    geom_point(shape = 21, color = "black", fill = "gray90", size = 4) +
    ggrepel::geom_text_repel(size = 6, color = "black", point.padding = unit(5, "lines")) +
    labs(x = "# Clusters", y = "Modularity") +
    GeneralTheme(18)

  b <- clustree::clustree(object,
                          prefix = "walktrap_knn",
                          node_size = 10,
                          node_text_colour = "black",
                          node_colour = "sc3_stability",
                          layout = "sugiyama",
                          edge_arrow = FALSE)
  b <- b + scale_color_continuous(low = "gray90", high = "#517AC9", name = "SC3 Stability") +
    ggraph::scale_edge_color_continuous(low = "gray90", high = "#517AC9", name = "Cell Count") +
    ggraph::scale_edge_alpha_continuous(name = "Cell Proportion")

  object$umap1 <- object@reductions$umap@cell.embeddings[, 1]
  object$umap2 <- object@reductions$umap@cell.embeddings[, 2]

  grid <- a + b + plot_layout(widths = c(0.5, 1))
  return(grid)

}

#' Define markers with `presto` with additional metrics
#'
#' @param seurat.obj Seurat object
#' @param group Column in `seurat.obj[[]]` to compare
#'
#' @return A data.frame of markers
#' @export
DefineMarkers <- function(seurat.obj, group) {
  markers <- presto::wilcoxauc(seurat.obj, group)
  n <- length(unique(markers$group))
  filtered_markers <- markers %>% dplyr::filter(padj <= 1E-3 & auc >= 0.5 & logFC >= log(1.1))
  StartFutureLapply()
  usethis::ui_info("Calculating gene specificity indices and max-to-second ratios")
  metrics <- future.apply::future_apply(filtered_markers, 1, function(x, markers_table) {
    tbl <- markers_table[markers_table$feature == x["feature"], ]
    max <- MaxN(tbl$avgExpr, 1)
    secondmax <- MaxN(tbl$avgExpr, 2)
    mtsr <- tbl$avgExpr[max]/tbl$avgExpr[secondmax]
    gsi <- tbl$avgExpr[tbl$group == x["group"]]/(1/n*sum(tbl$avgExpr[tbl$group != x["group"]]))
    return(data.frame(mtsr = mtsr, gsi = gsi))
  }, markers_table = markers)
  StopFutureLapply()
  metrics <- data.table::rbindlist(metrics)
  filtered_markers <- cbind(filtered_markers, metrics)
  return(filtered_markers)
}

MaxN <- function(x, n) {
  order(x, decreasing = TRUE)[n]
}

GetSTM <- function(x, idx) {
  x <- as.numeric(x[1:idx])
  max1 <- MaxN(x[1:idx], n = 1)
  max2 <- MaxN(x[1:idx], n = 2)
  diff <- x[max1]/x[max2]
  return(diff)
}

# Modified from Seurat
GroupSingletons <- function(ids, SNN, min.size = 9, clusters.to.merge, group.singletons = TRUE, verbose = TRUE) {

  # identify singletons
  singletons <- c()
  singletons <- names(x = which(x = table(ids) <= min.size))
  singletons <- intersect(x = unique(x = ids), singletons)

  if (!missing(clusters.to.merge)) {
    singletons <- c(singletons, as.character(clusters.to.merge))
  }

  if (!group.singletons) {
    ids[which(ids %in% singletons)] <- "singleton"
    return(ids)
  }

  # calculate connectivity of singletons to other clusters, add singleton
  # to cluster it is most connected to

  if (!is_empty(singletons)) {
    cluster_names <- as.character(x = unique(x = ids))
    cluster_names <- setdiff(x = cluster_names, y = singletons)
    connectivity <- vector(mode = "numeric", length = length(x = cluster_names))
    names(x = connectivity) <- cluster_names
    new.ids <- ids
    for (i in singletons) {
      i.cells <- names(which(ids == i))
      for (j in cluster_names) {
        j.cells <- names(which(ids == j))
        subSNN <- SNN[i.cells, j.cells]
        set.seed(1) # to match previous behavior, random seed being set in WhichCells
        if (is.object(x = subSNN)) {
          connectivity[j] <- sum(subSNN) / (nrow(x = subSNN) * ncol(x = subSNN))
        } else {
          connectivity[j] <- mean(x = subSNN)
        }
      }
      m <- max(connectivity, na.rm = T)
      mi <- which(x = connectivity == m, arr.ind = TRUE)
      closest_cluster <- sample(x = names(x = connectivity[mi]), 1)
      ids[i.cells] <- closest_cluster
    }
  }

  if (length(x = singletons) > 0 && verbose) {
    message(paste(
      length(x = singletons),
      "singletons identified.",
      length(x = unique(x = ids)),
      "final clusters."
    ))
  }

  return(ids)
}
