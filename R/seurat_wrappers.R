#' Process Seurat seurat.object through PCA
#'
#' @param seurat.obj Seurat object
#' @param n.var.features # of variable features to use
#' @param n.pcs # of PCs to use
#' @param remove.genes Genes to remove from variable features and PCA
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA
ReduceDims <- function(
  seurat.obj,
  n.var.features = 3000,
  n.pcs = 50,
  remove.genes = NULL,
  pca.name = "pca",
  ...
) {
  seurat.obj <- NormalizeData(object = seurat.obj)
  seurat.obj <- FindVariableFeatures(object = seurat.obj, selection.method = "vst", nfeatures = n.var.features)
  if (!missing(remove.genes) & !is.null(remove.genes)) {
    usethis::ui_done("Removed {length(remove.genes)} genes from variable features")
    VariableFeatures(seurat.obj) <- setdiff(VariableFeatures(seurat.obj), remove.genes)
  }
  seurat.obj <- ScaleData(object = seurat.obj, features = VariableFeatures(seurat.obj), block.size = 1000)
  usethis::ui_todo("Running PCA...")
  seurat.obj <- RunPCA(object = seurat.obj, npcs = n.pcs, features = VariableFeatures(seurat.obj),
                       verbose = FALSE, seed.use = 1, weight.by.var = TRUE, reduction.name = pca.name)
  return(seurat.obj)
}

#' Calculate heuristic metrics to determine optimal number of principal components for downstream analyses
#'
#' @param seurat.obj Seurat object
#' @param percent.stdev.range Percent of PCA standard deviation range to use for
#' @param set.pcs Number of PCs to use for maximum distance and derivative calculations
#' @param derivs.change.threshold Threshold used for second derivative changes
#'
#' @return Returns `object` with optimal PC information within the `object@misc` slot
#' @export
CalcPCHeuristics <- function(
  seurat.obj,
  reduction = "pca",
  store.name = "pca_metrics",
  percent.stdev.range = 0.05,
  set.pcs = NULL,
  force_tp = FALSE
) {

  middle_pts <- function(x) x[-1] - diff(x) / 2
  if (reduction == "harmony") {
    seurat.obj@reductions$harmony@stdev <- as.numeric(apply(seurat.obj@reductions$harmony@cell.embeddings, 2, sd))
  }

  object_stdev <- Seurat::Stdev(seurat.obj, reduction = reduction)
  stdev_range <- range(object_stdev)[2] - range(object_stdev)[1]
  cutoff <- min(object_stdev) + stdev_range * percent.stdev.range

  if (is.null(set.pcs)) {
    pcs_to_use <- max(which(object_stdev > cutoff))
  } else {
    pcs_to_use <- set.pcs
  }

  pcs_stdev <- object_stdev[1:pcs_to_use]
  pcs <- seq(1:length(pcs_stdev))

  slope <- (pcs_stdev[1] - pcs_stdev[length(pcs)])/(pcs[1]-pcs[length(pcs)])
  intercept <- pcs_stdev[1]-slope*pcs[1]
  b <- c(pcs[1], pcs_stdev[1])
  c <- c(pcs[length(pcs_stdev)], pcs_stdev[length(pcs_stdev)])
  dist <- vector()
  for (i in seq_along(1:length(pcs_stdev))) {
    a <- c(pcs[i], pcs_stdev[i])
    v1 <- b - c
    v2 <- a - b
    m <- cbind(v1,v2)
    dist[i] <- abs(det(m))/sqrt(sum(v1*v1))
  }

  derivs <- list()
  derivs$max_dist_pcs <- which.max(dist)
  derivs$percent_cutoff_pcs <- max(which(object_stdev > cutoff))
  derivs$slope <- slope
  derivs$intercept <- intercept

  if (ncol(seurat.obj) < 50000 | force_tp == TRUE) {
    tp <- intrinsicDimension::maxLikGlobalDimEst(seurat.obj@reductions$pca@cell.embeddings, k = 10)
    tp <- ceiling(tp$dim.est)
  } else {
    tp <- 101
  }
  derivs$tp_pc <- tp
  seurat.obj@misc[[store.name]] <- derivs
  return(seurat.obj)
}


Embed <- function(
  seurat.obj,
  reduction = "pca",
  heuristics.name = "pca_metrics",
  dims.use,
  knn.use = 20,
  store.graph = TRUE
) {
  dims_name <- glue::glue("{reduction}_dims")
  if (!missing(dims.use)) {
    seurat.obj@misc[[dims_name]] <- dims.use
  } else if (unique(seurat.obj@misc[[heuristics.name]]$tp_pc) == 101 | unique(seurat.obj@misc[[heuristics.name]]$tp_pc) == 0 | is.nan(unique(seurat.obj@misc[[heuristics.name]]$tp_pc))) {
    seurat.obj@misc[[dims_name]] <- unique(seurat.obj@misc[[heuristics.name]]$percent_cutoff_pcs)
  } else {
    seurat.obj@misc[[dims_name]] <- unique(seurat.obj@misc[[heuristics.name]]$tp_pc)
  }
  usethis::ui_info("\nUsing {seurat.obj@misc[[dims_name]]} dimensions for UMAP...")

  seurat.obj <- Seurat::RunUMAP(object = seurat.obj, reduction = reduction, reduction.name = glue::glue("{reduction}_umap"),
                                reduction.key = glue::glue("{substring(reduction, 1, 1)}UMAP_"),
                                dims = 1:seurat.obj@misc[[dims_name]], seed.use = 1)
  if (store.graph) {
    seurat.obj <- PrepareGraph(object = seurat.obj, reduction = reduction, dims = seurat.obj@misc[[dims_name]], knn = knn.use)
  }
  return(seurat.obj)
}

Walktrap <- function(
  seurat.obj,
  col.name = "pca_walktrap",
  graph.name = "pca_snn_20", ...
) {
  walktrap <- igraph::cluster_walktrap(seurat.obj@misc[[graph.name]], ...)
  ids <- walktrap$membership
  names(ids) <- Cells(seurat.obj)
  ids <- GroupSingletons(ids, seurat.obj@graphs[[graph.name]],
                         min.size = 9, group.singletons = TRUE, verbose = TRUE)
  seurat.obj <- AddMetaData(seurat.obj, ids, col.name = col.name)
  usethis::ui_done("Identified {length(unique(seurat.obj[[col.name, drop = TRUE]]))} walktrap clusters")
  return(seurat.obj)
}


IdentifyControlFeatures <- function(object, features,
                                    pool= NULL,
                                    nbin = 20,
                                    ctrl = 100,
                                    seed = NULL,
                                    assay = NULL) {
  # set seed
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }

  # extract assay data and find rowmeans
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay
  assay.data <- GetAssayData(object = object)
  rowmeans <- Matrix::rowMeans(assay.data, na.rm = TRUE)
  assay.data <- assay.data[rowmeans > 0, ]
  features.old <- features

  # check input features
  if (is.null(x = features)) {
    stop("Missing input feature list")
  }
  features <- lapply(
    features, function(x) {
      return(intersect(x = x, y = rownames(x = assay.data)))
    }
  )
  cluster.length <- length(x = features)
  if (!all(Seurat:::LengthCheck(values = features))) {
    warning(paste(
      'Could not find enough features in the object from the following feature lists:',
      paste(names(x = which(x = !Seurat:::LengthCheck(values = features)))),
      'Attempting to match case...'
    ))
    features <- lapply(
      X = features.old,
      FUN = CaseMatch,
      match = rownames(x = object)
    )
  }
  if (!all(Seurat:::LengthCheck(values = features))) {
    stop(paste(
      'The following feature lists do not have enough features present in the object:',
      paste(names(x = which(x = !Seurat:::LengthCheck(values = features)))),
      'exiting...'
    ))
  }

  pool <- rownames(x = assay.data)
  data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- ggplot2::cut_number(x = data.avg + stats::rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(x = sample(
          x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
          size = ctrl,
          replace = FALSE
        ))
      )
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  return(ctrl.use)
}

#' AddModuleScore from Seurat using scaled expression values
#'
#' See `?Seurat::AddModuleScore`
#'
#' @param object Seurat object
#' @param features Genes to use
#' @param pool Genes to pull from for controls
#' @param nbin Number of bins
#' @param ctrl Number of control genes
#' @param k Feature clusters
#' @param assay Assay to use
#' @param name Prefix for metadata columns
#' @param seed Seed
#' @param search Search for gene names
#' @param ...
#'
#' @return Seurat object
#' @export
#' @importFrom Seurat DefaultAssay `DefaultAssay<-` GetAssayData UpdateSymbolList CaseMatch
#' @importfrom
AddModuleScoreScaled <- function(
  object,
  features,
  pool = NULL,
  nbin = 20,
  ctrl = 100,
  k = FALSE,
  assay = NULL,
  name = "Cluster",
  seed = 1,
  ...
) {
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay
  assay.data <- GetAssayData(object = object)
  rowmeans <- Matrix::rowMeans(assay.data, na.rm = TRUE)
  assay.data <- assay.data[rowmeans > 0, ]
  features.old <- features
  if (k) {
    .NotYetUsed(arg = 'k')
    features <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      features[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == i))
    }
    cluster.length <- length(x = features)
  } else {
    if (is.null(x = features)) {
      stop("Missing input feature list")
    }
    features <- lapply(
      features, function(x) {
        return(intersect(x = x, y = rownames(x = assay.data)))
      }
    )
    cluster.length <- length(x = features)
  }
  if (!all(Seurat:::LengthCheck(values = features))) {
    warning(paste(
      'Could not find enough features in the object from the following feature lists:',
      paste(names(x = which(x = !Seurat:::LengthCheck(values = features)))),
      'Attempting to match case...'
    ))
    features <- lapply(
      X = features.old,
      FUN = CaseMatch,
      match = rownames(x = object)
    )
  }
  if (!all(Seurat:::LengthCheck(values = features))) {
    stop(paste(
      'The following feature lists do not have enough features present in the object:',
      paste(names(x = which(x = !Seurat:::LengthCheck(values = features)))),
      'exiting...'
    ))
  }
  pool <- rownames(x = assay.data)
  data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- ggplot2::cut_number(x = data.avg + stats::rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(x = sample(
          x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
          size = ctrl,
          replace = FALSE
        ))
      )
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(x = ctrl.use),
    ncol = ncol(x = object)
  )

  all_genes <- c(unique(unlist(ctrl.use)), unique(unlist(features)))
  usethis::ui_todo("Scaling data...")
  assay.data <- t(apply(assay.data[all_genes, ], 1, scale, center = FALSE))
  usethis::ui_done("Done")
  Clean()
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, ])
  }

  features.scores <- matrix(
    data = numeric(length = 1L),
    nrow = cluster.length,
    ncol = ncol(x = object)
  )
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- assay.data[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  assertthat::assert_that(!all(is.nan(features.scores)))
  assertthat::assert_that(!all(is.nan(ctrl.scores)))
  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = object)
  object[[colnames(x = features.scores.use)]] <- features.scores.use
  Seurat:::CheckGC()
  DefaultAssay(object = object) <- assay.old
  return(object)
}

FilterObject <- function(object,
                         nmads = 5,
                         variables = c("nCount_RNA", "nFeature_RNA", "percent_ribo", "percent_mito"),
                         batch = "batch",
                         percent_mito_cutoff = 20,
                         percent_ribo_cutoff = 50,
                         percent_hsp_cutoff = 10,
                         percent_soup_cutoff = 50,
                         include_soup = FALSE,
                         nCount_RNA_cutoff = 200,
                         nFeature_RNA_cutoff = 100,
                         remove = FALSE,
                         scrub = FALSE,
                         qc_column = "passing_qc") {

  threshold_drop <- c()
  threshold_drop <- c(threshold_drop, colnames(object)[object$percent_mito > percent_mito_cutoff])
  threshold_drop <- c(threshold_drop, colnames(object)[object$percent_ribo > percent_ribo_cutoff])
  threshold_drop <- c(threshold_drop, colnames(object)[object$percent_hsp > percent_hsp_cutoff])
  if (include_soup) {
    threshold_drop <- c(threshold_drop, colnames(object)[object$percent_soup > percent_soup_cutoff])
  }
  threshold_drop <- c(threshold_drop, colnames(object)[object$nCount_RNA < nCount_RNA_cutoff])
  threshold_drop <- c(threshold_drop, colnames(object)[object$nFeature_RNA < nFeature_RNA_cutoff])

  threshold_drop <- unique(threshold_drop)
  usethis::ui_info("{length(threshold_drop)} cells dropped by thresholds")
  keep_cells <- setdiff(Cells(object), threshold_drop)
  thresholded_object <- subset(object, cells = keep_cells)
  head(thresholded_object[[]])

  batch_use <- as.factor(thresholded_object$batch)
  mad_drop <- map(variables, ~ {
    metric_use <- thresholded_object[[.x, drop = TRUE]]
    dropped <- Cells(thresholded_object)[scater::isOutlier(metric = metric_use, nmads = nmads,
                                                           log = TRUE, type = "both", batch = batch_use)]
    return(dropped)
  })
  all_drop <- unique(c(unlist(mad_drop), threshold_drop))
  all_drop <- all_drop[!is.na(all_drop)]
  usethis::ui_info("{length(all_drop)} cells dropped by all filters")

  if (remove) {
    object <- object[, !(all_drop)]
  } else {
    object[[qc_column]] <- TRUE
    object[[qc_column]][all_drop, ] <- FALSE
  }

  if (scrub) {
    scrubbed <- Scrub(object[, setdiff(Cells(object), all_drop)], batch.var = "batch")
    scrubbed_scores <- scrubbed[[]] %>% select(scrublet_score, scrublet_label)
    object <- AddMetaData(object, scrubbed_scores)
    object$scrublet_score <- as.numeric(object$scrublet_score)

    summary(object$scrublet_score[object$scrublet_label == TRUE])
    summary(object$scrublet_score[object$scrublet_label == FALSE])
    object$passing_qc[object$scrublet_score > 0.05] <- FALSE
  }

  return(object)
}

RunAll <- function(seurat.obj, batch.var = "batch", nolist, harmony.dims = 30, ...) {

  seurat.obj <- ReduceDims(seurat.obj, n.var.features = 3000, n.pcs = 50, remove.genes = nolist)
  seurat.obj <- CalcPCHeuristics(seurat.obj, reduction = "pca", store.name = "pca_metrics", force_tp = TRUE)
  usethis::ui_info(seurat.obj@misc$pca_metrics$tp_pc)
  Clean()
  seurat.obj <- harmony::RunHarmony(seurat.obj, group.by.vars = "batch", dims.use = 1:harmony.dims,
                                    max.iter.harmony = 100, verbose = TRUE, assay.use = "RNA", plot_convergence = FALSE,
                                    theta = 1, sigma = 0.1, lambda = 1, max.iter.cluster = 30)
  Clean()
  seurat.obj <- CalcPCHeuristics(seurat.obj, reduction = "harmony", store.name = "harmony_metrics", force_tp = TRUE)
  usethis::ui_info(seurat.obj@misc$harmony_metrics$tp_pc)
  seurat.obj <- Embed(seurat.obj = seurat.obj, reduction = "harmony", heuristics.name = "harmony_metrics", ...)
  Clean()
  seurat.obj <- Walktrap(seurat.obj = seurat.obj, col.name = "harmony_walktrap", graph.name = "harmony_snn_20")
  Clean()
  StartFuture()
  seurat.obj <- AutoCluster(seurat.obj = seurat.obj, graph.name = "harmony_snn_20", col.name = "harmony_leiden", mod.percentile = 1)
  Clean()
  StopFuture()

  seurat.obj <- SetIdent(seurat.obj, value = "harmony_leiden")
  seurat.obj <- BuildClusterTree(seurat.obj, reorder = TRUE, reorder.numeric = TRUE)
  seurat.obj[["harmony_leiden"]] <- Idents(seurat.obj)
  return(seurat.obj)

}


AddBatchFreq <- function(object, batch = "batch") {
  sample_df <- as.data.frame(table(object[[batch]]))
  colnames(sample_df) <- c("batch", "freq")
  merged <- merge(object[[]], sample_df, by = "batch")
  rownames(merged) <- merged$cell_barcode
  merged <- merged %>% select(freq)
  object <- AddMetaData(object, merged, col.name = "batch_frequency")
  return(object)
}


#' Run Scrublet
#'
#' @param seurat.obj Seurat object
#' @param batch.var Batch variable
#'
#' @return Seurat object
#' @export
Scrub <- function(seurat.obj, batch.var = "batch") {
  scrublet <- reticulate::import("scrublet")
  seurat.obj$scrublet_score <- "NA"
  seurat.obj$scrublet_label <- "NA"
  sample_df <- as.data.frame(table(seurat.obj[[batch.var]]))
  colnames(sample_df) <- c("batch", "freq")
  sample_df$batch <- as.character(sample_df$batch)
  sample_df$freq <- as.numeric(sample_df$freq)

  for (i in seq_along(1:length(sample_df$batch))) {
    freq <- as.numeric(sample_df$freq[i])
    cells <- (seurat.obj[[batch.var, drop = TRUE]] == sample_df$batch[i])
    if (freq < 100) {
      message(glue(">> Only {freq} cells, skipping doublet prediction"))
      seurat.obj[["scrublet_score"]][cells, ] <- NA
      seurat.obj[["scrublet_label"]][cells, ] <- NA
    } else {
      matrix <- as.matrix(GetAssayData(seurat.obj, slot = "counts")[, cells])
      scrublet_object <- scrublet$Scrublet(t(matrix), expected_doublet_rate = 4.6e-06*freq)
      message(glue(">> Scrublet object created for iteration {i}/{length(sample_df$batch)}"))
      scores <- scrublet_object$scrub_doublets(min_counts = 3, min_cells = 3,
                                               min_gene_variability_pctl = 85, n_prin_comps = as.integer(30), verbose = TRUE)
      message(glue(">> Identified {sum(as.vector(scores[[2]]))}/{length(scores[[2]])} cells as doublets"))
      seurat.obj[["scrublet_score"]][cells, ] <- scores[[1]]
      seurat.obj[["scrublet_label"]][cells, ] <- scores[[2]]
    }
  }
  return(seurat.obj)
}

#' Add percentage expression families
#' @param object Seurat object to add metadata to
#' @return Seurat object with metadata columns added
AddExprMeta <- function(object) {

  object[["percent_mito"]] <- PercentageFeatureSet(object, pattern = "^MT-|^MTRNR|^MTERF|^MTFR")

  # https://www.genenames.org/data/genegroup/#!/group/1054
  ribo <- intersect(rownames(object), c("MRPL1", "MRPL2", "MRPL3", "MRPL4", "MRPL10", "MRPL11", "MRPL12",
                                        "MRPL13", "MRPL14", "MRPL15", "MRPL16", "MRPL17", "MRPL18", "MRPL19",
                                        "MRPL20", "MRPL21", "MRPL22", "MRPL23", "MRPL24", "MRPL27", "MRPL28",
                                        "MRPL30", "MRPL32", "MRPL33", "MRPL34", "MRPL35", "MRPL36", "MRPL37",
                                        "MRPL38", "MRPL39", "MRPL40", "MRPL41", "MRPL42", "MRPL43", "MRPL44",
                                        "MRPL45", "MRPL46", "MRPL47", "MRPL48", "MRPL49", "MRPL50", "MRPL51",
                                        "MRPL52", "MRPL53", "MRPL54", "MRPL55", "MRPL57", "MRPL58", "RPLP0",
                                        "RPLP1", "RPLP2", "RPL3", "RPL3L", "RPL4", "RPL5", "RPL6", "RPL7",
                                        "RPL7A", "RPL7L1", "RPL8", "RPL9", "RPL10", "RPL10A", "RPL10L",
                                        "RPL11", "RPL12", "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18A",
                                        "RPL19", "RPL21", "RPL22", "RPL23", "RPL23A", "RPL24", "RPL26",
                                        "RPL26L1", "RPL27", "RPL27A", "RPL28", "RPL29", "RPL30", "RPL31",
                                        "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", "RPL36A", "RPL36AL",
                                        "RPL37", "RPL37A", "RPL38", "RPL39", "RPL39L", "UBA52", "RPL41",
                                        "MRPL1", "MRPL2", "MRPL3", "MRPL4", "MRPL9", "MRPL10", "MRPL11",
                                        "MRPL12", "MRPL13", "MRPL14", "MRPL15", "MRPL16", "MRPL17", "MRPL18",
                                        "MRPL19", "MRPL20", "MRPL21", "MRPL22", "MRPL23", "MRPL24", "MRPL27",
                                        "MRPL28", "MRPL30", "MRPL32", "MRPL33", "MRPL34", "MRPL35", "MRPL36",
                                        "MRPL37", "MRPL38", "MRPL39", "MRPL40", "MRPL41", "MRPL42", "MRPL43",
                                        "MRPL44", "MRPL45", "MRPL46", "MRPL47", "MRPL48", "MRPL49", "MRPL50",
                                        "MRPL51", "MRPL52", "MRPL53", "MRPL54", "MRPL55", "MRPL57", "MRPS2",
                                        "MRPS5", "MRPS6", "MRPS7", "MRPS9", "MRPS10", "MRPS11", "MRPS12",
                                        "MRPS14", "MRPS15", "MRPS16", "MRPS17", "MRPS18A", "MRPS18B",
                                        "MRPS18C", "MRPS21", "MRPS22", "MRPS23", "MRPS24", "MRPS25",
                                        "MRPS26", "MRPS27", "MRPS28", "DAP3", "MRPS30", "MRPS31", "MRPS33",
                                        "MRPS34", "MRPS35", "MRPS36", "MRPS2", "MRPS10", "MRPS11", "MRPS12",
                                        "MRPS14", "MRPS15", "MRPS16", "MRPS17", "MRPS18A", "MRPS18B",
                                        "MRPS18C", "MRPS21", "MRPS22", "MRPS23", "MRPS24", "MRPS25",
                                        "MRPS26", "MRPS27", "MRPS28", "MRPS30", "MRPS31", "MRPS33", "MRPS34",
                                        "MRPS35", "MRPS36", "RPS2", "RPS3", "RPS3A", "RPS4X", "RPS4Y1",
                                        "RPS4Y2", "RPS5", "RPS6", "RPS7", "RPS8", "RPS9", "RPS10", "RPS11",
                                        "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS17",
                                        "RPS18", "RPS19", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25",
                                        "RPS26", "RPS27", "RPS27A", "RPS27L", "RPS28", "RPS29", "FAU"))
  object[["percent_ribo"]] <- PercentageFeatureSet(object, features = ribo)

  # https://www.genenames.org/data/genegroup/#!/group/582
  hsp <- intersect(rownames(object), c("BBS10", "BBS12", "TCP1", "CCT2", "CCT3", "CCT4", "CCT5", "CCT6A",
                                       "CCT6B", "CCT7", "CCT8", "CLPB", "HSPD1", "HSPE1", "MKKS", "DNAJA1",
                                       "DNAJA2", "DNAJA3", "DNAJA4", "DNAJB1", "DNAJB2", "DNAJB3", "DNAJB4",
                                       "DNAJB5", "DNAJB6", "DNAJB7", "DNAJB8", "DNAJB9", "DNAJB11",
                                       "DNAJB12", "DNAJB13", "DNAJB14", "DNAJC1", "DNAJC2", "DNAJC3",
                                       "DNAJC4", "DNAJC5", "DNAJC5B", "DNAJC5G", "DNAJC6", "DNAJC7",
                                       "DNAJC8", "DNAJC9", "DNAJC10", "DNAJC11", "DNAJC12", "DNAJC13",
                                       "DNAJC14", "DNAJC15", "DNAJC16", "DNAJC17", "DNAJC18", "DNAJC19",
                                       "HSCB", "DNAJC21", "DNAJC22", "SEC63", "DNAJC24", "DNAJC25",
                                       "GAK", "DNAJC27", "DNAJC28", "SACS", "DNAJC30", "HSPA1A", "HSPA1B",
                                       "HSPA1L", "HSPA2", "HSPA4", "HSPA4L", "HSPA5", "HSPA6", "HSPA7",
                                       "HSPA8", "HSPA9", "HSPA12A", "HSPA12B", "HSPA13", "HSPA14", "HSPH1",
                                       "HYOU1", "HSP90AA1", "HSP90AA3P", "HSP90AB1", "HSP90B1", "TRAP1",
                                       "HSPB1", "HSPB2", "HSPB3", "CRYAA", "CRYAB", "HSPB6", "HSPB7",
                                       "HSPB8", "HSPB9", "ODF1", "HSPB11"))

  object[["percent_hsp"]] <- PercentageFeatureSet(object, features = hsp)

  hb <- intersect(rownames(object),
                  c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"))
  object[["percent_hgb"]] <- PercentageFeatureSet(object, features = hb)

  return(object)
}

FilterGenes <- function(object, assay = "RNA") {
  genes_threshold <- ifelse(nrow(object) > 10000, 10, 5)
  genes_filter <- Matrix::rowSums(object@assays[[assay]]@counts > 0)
  keep_genes <- genes_filter > genes_threshold
  length(keep_genes)
  object@assays[[assay]]@counts <- object@assays[[assay]]@counts[keep_genes, ]
  object@assays[[assay]]@data <- object@assays[[assay]]@data[keep_genes, ]
  return(object)
}

RemoveLowFreqBatches <- function(object, threshold = 500) {
  object$batch <- as.character(object$batch)
  if (any(table(object$batch) < threshold)) {
    ui_info("Removing low frequency batches")
    passing_batches <- names(table(object$batch))[table(object$batch) >= threshold]
    object <- object[, object$batch %in% passing_batches]
  } else {
    ui_info("No low frequency batches")
  }
  return(object)
}

LabelCyclingCells <- function(object, grouping.var = "seurat_clusters") {

  object@meta.data[, grep("s_score", colnames(object@meta.data))] <- NULL
  object@meta.data[, grep("g2m_score", colnames(object@meta.data))] <- NULL
  object@meta.data[, grep("phase", colnames(object@meta.data))] <- NULL

  # calculate cell cycle scores
  object <- Seurat::CellCycleScoring(object,
                                     s.features = cc.genes.updated.2019$s.genes,
                                     g2m.features = cc.genes.updated.2019$g2m.genes,
                                     set.ident = FALSE, search = FALSE)
  head(object[[]])
  colnames(object@meta.data)[grep("S.Score", colnames(object@meta.data))] <- "s_score"
  colnames(object@meta.data)[grep("G2M.Score", colnames(object@meta.data))] <- "g2m_score"
  colnames(object@meta.data)[grep("Phase", colnames(object@meta.data))] <- "phase"
  object$cc_diff <- object$s_score - object$g2m_score

  # summarize scores based on clusters
  summary <- object@meta.data %>%
    group_by(!!sym(grouping.var)) %>%
    summarize(median = median(cc_diff), mean = mean(cc_diff), sd = sd(cc_diff))

  # calculate outlier cell cutoffs
  cc_diff_norm_mean <- mean(summary$mean[  !(summary$mean > (median(summary$mean) + 2*mad(summary$mean))) | ! (summary$mean <(median(summary$mean) - 2*mad(summary$mean)))])
  cc_diff_norm_sd <- mean(summary$sd[  !(summary$mean > (median(summary$mean) + 2*mad(summary$mean))) | ! (summary$mean <(median(summary$mean) - 2*mad(summary$mean)))])
  upper <- cc_diff_norm_mean + 2*cc_diff_norm_sd
  lower <- cc_diff_norm_mean - 2*cc_diff_norm_sd

  # assign labels
  object$cc_diff_label <- ifelse(object$cc_diff > upper | object$cc_diff < lower, "Cycling", "Non-cycling")

  return(object)
}
