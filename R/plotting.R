#' Save png, pdf and plot object
#'
#' @param plot `ggplot2` object to save
#' @param filename Filename
#' @param root Directories to add to /plots or /data
#' @param h Height
#' @param w Width
#' @param s Scale
#'
#' @return NULL
#' @export
SavePlot <- function(
  plot,
  filename,
  root = "publication",
  h = 6,
  w = 6,
  s = 1,
  save.data = TRUE
) {
  ggplot2::ggsave(plot = plot, filename = glue::glue("plots/{root}/{filename}.png"),
                  scale = s, width = w, height = h, units = "in", dpi = 300)
  ggplot2::ggsave(plot = plot, filename = glue::glue("plots/{root}/{filename}.pdf"),
                  scale = s, width = w, height = h, units = "in", dpi = 300)
  if (save.data) {
    saveRDS(plot, file = glue::glue("data/{root}/plots/{filename}.rds"))
  }
  usethis::ui_done("Saved")
}

#' RemoveAxes
#'
#' Modified from Seurat::NoAxes()
#'
#' @return
#' @export
RemoveAxes <- function (..., keep.text = FALSE, keep.ticks = FALSE)
{
  blank <- element_blank()
  no_axes_theme <- theme(axis.line.x = blank,
                         axis.line.y = blank,
                         validate = TRUE, ...)
  if (!keep.text) {
    no_axes_theme <- no_axes_theme + theme(
      axis.text.x = blank,
      axis.text.y = blank,
      validate = TRUE, ...)
  }

  if (!keep.ticks) {
    no_axes_theme <- no_axes_theme + theme(
      axis.ticks.x = blank,
      axis.ticks.y = blank,
      validate = TRUE, ...)
  }

  return(no_axes_theme)
}

#' Remove backgrounds
#'
#' @param outline Keep plot outline
#' @param ...
#'
#' @return
#' @export
RemoveBackgrounds <- function(outline = FALSE, ...)
{
  if (outline) {
    no_bg_theme <- theme(panel.background = element_rect(fill = "transparent", color = "black", size = 1),
                         plot.background = element_rect(fill = "transparent", color = "transparent", size = 0),
                         legend.background = element_rect(fill = "transparent", color = "transparent", size = 0),
                         legend.box.background = element_rect(fill = "transparent", color = "transparent", size = 0),
                         panel.grid = element_blank(),
                         axis.line = element_blank(),
                         validate = TRUE, ...)
  } else {
    no_bg_theme <- theme(panel.background = element_rect(fill = "transparent", color = "transparent", size = 0),
                         plot.background = element_rect(fill = "transparent", color = "transparent", size = 0),
                         legend.background = element_rect(fill = "transparent", color = "transparent", size = 0),
                         legend.box.background = element_rect(fill = "transparent", color = "transparent", size = 0),
                         panel.grid = element_blank(),
                         validate = TRUE, ...)
  }

  return(no_bg_theme)
}


#' Space axis titles away from plot area
#'
#' @param scale Increase spacing
#' @param ...
#'
#' @return
#' @export
SpaceAxisTitles <- function(scale = 1, ...) {
  theme (
    axis.title.x = element_text(face = "plain", margin = margin(12*scale, 0, 0, 0)),
    axis.title.y = element_text(face = "plain", margin = margin(0, 12*scale, 0, 0)),
    validate = TRUE
  )
}

#' General plotting theme
#'
#' @param base_size
#' @param ...
#'
#' @return
#' @export
GeneralTheme <- function(base_size, ...) {
  theme_classic(base_size = base_size, ...) +
    ggeasy::easy_all_text_color("black") +
    theme(
      axis.line = element_blank(),
      plot.title = element_text(size =  base_size, color = "black", face = "bold", margin = margin(0,0,4,0)),
      plot.subtitle = element_text(size = base_size - 2, color = "black", margin = margin(0,0,4,0)),
      panel.background = element_rect(fill = "transparent", color = "black", size = 1),
      plot.background = element_rect(fill = "transparent", color = "transparent", size = 0),
      panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
      plot.caption = element_text(hjust = 0, color = "gray40", margin = margin(12)),
      legend.title = element_text(size = base_size, face = "plain"),
      legend.text = element_text(size = base_size - 2),
      legend.background = element_rect(fill = "transparent", color = "transparent"),
      legend.box.background = element_rect(fill = "transparent", color = "transparent"),
      legend.position = "right",
      legend.justification = "top",
      legend.key.size = unit(1, "line"),
      validate = TRUE
    )
}

PlotLabelOverlap <- function(object, xvar, yvar) {

  res_props <- object[[]] %>%
    dplyr::group_by(!!sym(xvar), !!sym(yvar)) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(prop = n / sum(n))
  res_props <- res_props %>% dplyr::select(-n) %>% spread(!!sym(yvar), prop)
  res_props[is.na(res_props)] <- 0
  rows <- res_props[, 1, drop = TRUE]
  res_props <- res_props[, -1]
  res_props <- data.matrix(res_props)
  rownames(res_props) <- rows

  cols <- rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")[c(1, 4, 7)])
  col_fun <- circlize::colorRamp2(c(0, 1), cols[c(2,3)])
  ht <- ComplexHeatmap::Heatmap(
    res_props,
    col = col_fun,
    row_names_side = "right",
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    show_column_dend = FALSE,
    show_row_dend = FALSE,
    row_names_gp = grid::gpar(fontface = "bold", fontsize = 12),
    row_names_centered = FALSE,
    column_names_gp = grid::gpar(fontface = "bold", fontsize = 12),
    column_names_centered = FALSE,
    column_names_rot = 45,
    border_gp = grid::gpar(col = "black", lwd = 2),
    border = TRUE,
    rect_gp = grid::gpar(col = "black", lwd = 1),
    row_dend_width = unit(32, "pt"),
    heatmap_legend_param = list(title = "Fractional\nabundance\n",
                                border = "black", legend_height = unit(128, "pt"),
                                legend_width = unit(128, "pt"), direction = "vertical", legend_border = grid::gpar(col = "black", lwd = 2),
                                labels_gp = grid::gpar(fontsize = 12), title_gp = grid::gpar(fontsize = 12, fontface = "plain")))
  return(ht)
}



#' Modified version of `Seurat::Dotplot` with outlined points
#'
#' See `?Seurat::DotPlot`
#' @param object Object
#' @param assay Assay
#' @param features Features to plot
#' @param cols cols
#' @param col.min col.min
#' @param col.max col.max
#' @param dot.min dot.min
#' @param dot.scale dot.scale
#' @param group.by group.by
#' @param split.by split.by
#' @param scale.by scale.by
#' @param scale.min scale.min
#' @param scale.max scale.max
#'
#' @return ggplot2 object
#' @export
OutlinedDotPlot <- function (object, assay = NULL, features, cols = c("lightgrey",
                                                                      "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
                             group.by = NULL, split.by = NULL, scale.by = "radius", scale.min = NA,
                             scale.max = NA)
{
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  scale.func <- switch(EXPR = scale.by, size = scale_size,
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  data.features <- FetchData(object = object, vars = features)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)
  }
  else {
    object[[group.by, drop = TRUE]]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]]
    if (length(x = unique(x = splits)) > length(x = cols)) {
      stop("Not enought colors for the number of groups")
    }
    cols <- cols[1:length(x = unique(x = splits))]
    names(x = cols) <- unique(x = splits)
    data.features$id <- paste(data.features$id, splits,
                              sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)),
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident,
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = Seurat:::PercentAbove,
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot ==
                                                     x, "avg.exp"]
                             data.use <- scale(x = data.use)
                             data.use <- MinMax(data = data.use, min = col.min,
                                                max = col.max)
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (!is.null(x = split.by)) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled,
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot,
                                    levels = rev(x = features))
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (!is.null(x = split.by)) {
    splits.use <- vapply(X = strsplit(x = as.character(x = data.plot$id),
                                      split = "_"), FUN = "[[", FUN.VALUE = character(length = 1L),
                         2)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled",
                     no = "colors")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot",
                                                        y = "id")) +
    geom_point(mapping = aes_string(size = "pct.exp", fill = color.by), shape = 21, alpha = 0.9) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    guides(size = guide_legend(title = "Percent\nExpressed")) +
    labs(x = "Features", y = ifelse(test = is.null(x = split.by),
                                    yes = "Identity", no = "Split Identity")) + cowplot::theme_cowplot()
  if (!is.null(x = split.by)) {
    plot <- plot + scale_fill_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + scale_fill_distiller(palette = cols)
  }
  else {
    plot <- plot + scale_fill_gradient(low = cols[1], high = cols[2])
  }
  if (is.null(x = split.by)) {
    plot <- plot + guides(fill = guide_colorbar(title = "Average\nExpression"))
  }
  return(plot)
}

PlotOverview <- function(seurat.obj, cluster.name, reduction.name, overview.name) {

  seurat.obj <- SetIdent(seurat.obj, value = cluster.name)
  seurat.obj <- BuildClusterTree(seurat.obj, reorder = TRUE, reorder.numeric = TRUE)
  seurat.obj@meta.data[, cluster.name] <- Idents(seurat.obj)

  a <- DimPlot(seurat.obj, reduction = reduction.name, group.by = "animal_array_label",
               label = TRUE, repel = TRUE, shuffle = TRUE, label.size = 4) +
    GeneralTheme(14) +
    labs(x = "UMAP1", y = "UMAP2", title = "") +
    scale_color_manual(values = seurat.obj[[]] %>% select(animal_array_label, animal_array_label_cols) %>% arrange(animal_array_label) %>% distinct() %>% pull(animal_array_label_cols)) +
    NoLegend()
  a$layers[[2]]$aes_params$fontface <- "bold"
  plot(a)
  ggsave(plot = a, filename = glue::glue("plots/wk10s2_{overview.name}_umap_batch.png"), dpi = 300, height = 8, width = 8, units = "in")

  b <- DimPlot(seurat.obj, reduction = reduction.name, group.by = cluster.name,
               label = TRUE, repel = TRUE, label.box = TRUE, label.color = "black", label.size = 4) +
    GeneralTheme(14) +
    labs(x = "UMAP1", y = "UMAP2", title = "") +
    NoLegend()
  b$layers[[2]]$aes_params$fontface <- "bold"
  plot(b)
  ggsave(plot = b, filename = glue::glue("plots/wk10s2_{overview.name}_umap_leiden.png"), dpi = 300, height = 8, width = 8, units = "in")

  # define cluster order
  tree <- as.hclust(seurat.obj@tools$BuildClusterTree)
  ddtree <- as.dendrogram(tree)
  cluster_tree <- ggtree::ggtree(ddtree) + scale_x_reverse() + coord_flip()

  # reorder markers
  leiden_markers <- presto::wilcoxauc(seurat.obj, group_by = cluster.name, seurat_assay = "RNA", assay = "data")
  top_leiden_markers <- leiden_markers %>% group_by(group) %>% top_n(2, auc)
  top_leiden_markers$group <- factor(top_leiden_markers$group, levels = tree$labels[tree$order])
  markers_to_plot <- top_leiden_markers %>%
    arrange(group) %>%
    pull(feature) %>%
    unique()

  dotplot <- ggplot(leiden_markers %>% filter(feature %in% markers_to_plot & pct_in > 1) %>%
                      mutate(group = factor(group, levels = tree$labels[tree$order]),
                             feature = factor(feature, levels = markers_to_plot)),
                    aes(x = group, y = feature, fill = auc, size = pct_in)) +
    geom_point(shape = 21, color = "black") +
    scale_size_continuous(range = c(0.1, 5), name = "% Expr.") +
    scale_fill_gradient2(low = "#3B4CC0", mid = "#DDDDDD", high = "#B40426", name = "auROC", midpoint = 0.5) +
    scale_y_discrete(position = "right") +
    labs(x = "Cluster", y = "") +
    GeneralTheme(base_size = 14) +
    theme(
      axis.line  = element_blank(),
      axis.text.x = element_text(size = 14, angle = 0, vjust = 0.5, hjust = 0.5),
      axis.ticks = element_blank(),
      panel.grid.major = element_line(color = "gray90", size = 0.5)
    )

  c <- cluster_tree + plot_spacer() + dotplot +
    plot_layout(nrow = 3, ncol = 1, heights = c(0.25, -0.05, 1))

  ggsave(plot = c, filename = glue::glue("plots/wk10s2_{overview.name}_dotplot_leiden.png"), dpi = 300, height = 8, width = 8, units = "in")

  return(seurat.obj)

}


ClusteredAUROCDotPlot <- function(
  object,
  ident
) {

  # define cluster order
  Idents(object) <- ident
  object <- BuildClusterTree(object, reorder = TRUE, reorder.numeric = FALSE)
  object[[ident]] <- Idents(object)
  tree <- as.hclust(object@tools$BuildClusterTree)
  ddtree <- as.dendrogram(tree)
  cluster_tree <- ggtree::ggtree(ddtree) + scale_x_reverse() + coord_flip()

  # reorder markers
  leiden_markers <- presto::wilcoxauc(object, group_by = ident,
                                      seurat_assay = "RNA", assay = "data")
  top_leiden_markers <- leiden_markers %>% group_by(group) %>% top_n(2, auc)
  top_leiden_markers$group <- factor(top_leiden_markers$group, levels = tree$labels[tree$order])
  markers_to_plot <- top_leiden_markers %>%
    arrange(group) %>%
    pull(feature) %>%
    unique()

  dotplot <- ggplot(leiden_markers %>% filter(feature %in% markers_to_plot & pct_in > 1) %>%
                      mutate(group = factor(group, levels = tree$labels[tree$order]),
                             feature = factor(feature, levels = markers_to_plot)),
                    aes(x = group, y = feature, fill = auc, size = pct_in)) +
    geom_point(shape = 21, color = "black") +
    scale_size_continuous(range = c(0.1, 5), name = "% Expr.") +
    scale_fill_gradient2(low = "#3B4CC0", mid = "whitesmoke", high = "#B40426", name = "AUROC", midpoint = 0.5) +
    scale_y_discrete(position = "right") +
    labs(x = "Cluster", y = "") +
    GeneralTheme(base_size = 14) +
    SpaceAxisTitles() +
    RemoveBackgrounds(outline = TRUE) +
    theme(
      axis.line  = element_blank(),
      axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
      axis.ticks = element_blank(),
      panel.grid.major = element_line(color = "gray90", size = 0.5),
      plot.margin = margin(0, 0, 0, 0.5, unit = "in")
    )
  dotplot
  f <- cluster_tree + plot_spacer() + dotplot +
    plot_layout(nrow = 3, ncol = 1, heights = c(0.25, -0.05, 1))
  return(f)
}
