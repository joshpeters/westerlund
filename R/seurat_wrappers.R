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
  n.var.features = 2000,
  n.pcs = 50,
  remove.genes,
  ...
) {
  seurat.obj <- NormalizeData(object = seurat.obj)
  seurat.obj <- FindVariableFeatures(object = seurat.obj, selection.method = "vst", nfeatures = n.var.features)
  if (!missing(remove.genes) & !is.null(remove.genes)) {
    usethis::ui_done("Removed {length(remove.genes)} genes from variable features")
    VariableFeatures(seurat.obj) <- setdiff(VariableFeatures(seurat.obj), remove.genes)
  }
  seurat.obj <- ScaleData(object = seurat.obj, features = VariableFeatures(seurat.obj), block.size = 1000)
  seurat.obj <- RunPCA(object = seurat.obj, npcs = n.pcs, features = VariableFeatures(seurat.obj), verbose = FALSE, seed.use = 1, weight.by.var = TRUE)
  seurat.obj <- CalcPCHeuristics(seurat.obj = seurat.obj, force_tp = TRUE, ...)
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
  set.pcs = NA,
  force_tp = FALSE
) {

  middle_pts <- function(x) x[-1] - diff(x) / 2
  if (reduction == "harmony") {
    seurat.obj@reductions$harmony@stdev <- as.numeric(apply(seurat.obj@reductions$harmony@cell.embeddings, 2, sd))
  }

  object_stdev <- Seurat::Stdev(seurat.obj, reduction = reduction)
  stdev_range <- range(object_stdev)[2] - range(object_stdev)[1]
  cutoff <- min(object_stdev) + stdev_range * percent.stdev.range

  if (is.na(set.pcs)) {
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
  search = FALSE,
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
        missing.features <- setdiff(x = x, y = rownames(x = assay.data))
        if (length(x = missing.features) > 0) {
          warning(
            "The following features are not present in the object: ",
            paste(missing.features, collapse = ", "),
            ifelse(
              test = search,
              yes = ", attempting to find updated synonyms",
              no = ", not searching for symbol synonyms"
            ),
            call. = FALSE,
            immediate. = TRUE
          )
          if (search) {
            tryCatch(
              expr = {
                updated.features <- UpdateSymbolList(symbols = missing.features, ...)
                names(x = updated.features) <- missing.features
                for (miss in names(x = updated.features)) {
                  index <- which(x == miss)
                  x[index] <- updated.features[miss]
                }
              },
              error = function(...) {
                warning(
                  "Could not reach HGNC's gene names database",
                  call. = FALSE,
                  immediate. = TRUE
                )
              }
            )
            missing.features <- setdiff(x = x, y = rownames(x = assay.data))
            if (length(x = missing.features) > 0) {
              warning(
                "The following features are still not present in the object: ",
                paste(missing.features, collapse = ", "),
                call. = FALSE,
                immediate. = TRUE
              )
            }
          }
        }
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
