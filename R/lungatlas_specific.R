#' Preprocess object
#'
#' @param object
#' @param dims.use
#'
#' @return
#' @export
Preprocess <- function(
  object,
  dims.use = 20
  ) {

  object <- FilterGenes(object)
  object <- RemoveLowFreqBatches(object)
  object@meta.data[grep("percent", colnames(object@meta.data))] <- NULL
  object <- AddExprMeta(object)
  object <- FilterObject(object, nmads = 3, variables = c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo"),
                         batch = "batch", percent_mito_cutoff = 20, percent_ribo_cutoff = 50, percent_hsp_cutoff = 5,
                         nCount_RNA_cutoff = 200, nFeature_RNA_cutoff = 100, remove = FALSE, scrub = FALSE, qc_column = "passing_qc")
  assertthat::assert_that(all.equal(colnames(object), rownames(object[[]])))
  keep_cells <- Seurat::WhichCells(object, expression = passing_qc == TRUE)
  ui_info("Number of cells pre-filtering: {ncol(object)}")
  ui_info("Number of cells post-filtering: {length(keep_cells)}")
  object <- object[, keep_cells]

  object <- ReduceDims(seurat.obj = object, n.var.features = 3000, n.pcs = 50, remove.genes = NULL)
  object <- harmony::RunHarmony(object, group.by.vars = "batch", dims.use = 1:30,
                                max.iter.harmony = 100, verbose = TRUE, assay.use = "RNA", plot_convergence = FALSE,
                                theta = 1, sigma = 0.1, lambda = 1, max.iter.cluster = 30)
  object <- CalcPCHeuristics(object, reduction = "harmony", store.name = "harmony_metrics", force_tp = TRUE)
  object <- Embed(seurat.obj = object, reduction = "harmony", heuristics.name = "harmony_metrics",
                  dims.use = dims.use, knn.use = 20)
  object <- AutoCluster(seurat.obj = object, col.name = "harmony_leiden",
                        graph.name = "harmony_snn_20", mod.percentile = 0.99)
  usethis::ui_done("Saving object...")
  saveRDS(object, file = glue("data/{set}/{set}_preprocessed_object.rds"))
  return(object)
}

#' Transfer reference dataset labels
#'
#' @param object
#'
#' @return
#' @export
TransferLabels <- function(
  object
) {

  reference <- readRDS("data/reference/trkr20_reference.rds")
  mapped <- symphony::mapQuery(
    exp_query = Seurat::GetAssayData(object, slot = "data", assay = "RNA"),
    metadata_query = object@meta.data,
    vars = c("batch"),
    ref_obj = reference,
    do_normalize = FALSE,
    verbose = TRUE,
    do_umap = FALSE
  )
  mapped <-  symphony::knnPredict(mapped, reference, reference$meta_data$trkr20_type, k = 5)
  object$trkr20_type <- as.character(mapped$meta_data$cell_type_pred_knn)

  type_df <- read.csv("data/reference/trkr20_cellhierarchy.csv")
  meta_add <- merge(object[[]], type_df[, 1:4] %>% distinct(),  by = "trkr20_type", sort = FALSE, all.x = TRUE, all.y = FALSE)
  meta_add <- meta_add %>% select(trkr20_lineage, trkr20_class, trkr20_type, trkr20_ismnp, cell_barcode)
  rownames(meta_add) <- meta_add$cell_barcode
  meta_add$cell_barcode <- NULL
  object <- Seurat::AddMetaData(object, meta_add)
  DimPlot(object, reduction = "harmony_umap", group.by = "trkr20_type", label = TRUE, repel = TRUE) + NoLegend()

  object <- LabelCyclingCells(object, grouping.var = "harmony_leiden")

  trkr20_degs <- readRDS("data/reference/trkr20_DEGs.rds")
  trkr20_degs <- lapply(trkr20_degs, function(x) intersect(x, rownames(object)))
  object@meta.data[, grep("_score", colnames(object@meta.data))] <- NULL
  object <- Seurat::AddModuleScore(object, features = trkr20_degs, nbin = 20, ctrl = 30,
                           seed = 1, name = "trkr20_score", search = FALSE)
  colnames(object@meta.data)[grep("trkr20_score", colnames(object@meta.data))] <- paste0(names(trkr20_degs), "_score")

  object$mnp_diff_score <- object$mnp_score - object$nonmnp_score

  reference <- readRDS("data/reference/azimuth_reference.rds")
  mapped <- symphony::mapQuery(
    exp_query = GetAssayData(object, slot = "data", assay = "RNA"),
    metadata_query = object@meta.data,
    vars = c("batch"),
    ref_obj = reference,
    do_normalize = FALSE,
    verbose = TRUE,
    do_umap = FALSE
  )
  mapped <-  symphony::knnPredict(mapped, reference, reference$meta_data$celltype.l2, k = 5)
  object$azimuth_celltype_l2 <- as.character(mapped$meta_data$cell_type_pred_knn)
  return(object)
}

LabelCellTypes <- function(
  object
){

  mnp_scores <- object[[]] %>% select(harmony_leiden, mnp_diff_score)

  # fit mixture model
  mm <- mixtools::normalmixEM(mnp_scores$mnp_diff_score, k = 2)
  x_data <- with(mm, seq(min(x), max(x), len = 1000))
  pars <- with(mm, data.frame(comp = colnames(posterior), mu, sigma, lambda))
  em_df <- data.frame(x = rep(x_data, each = nrow(pars)), pars)
  em_df$y <- em_df$lambda * dnorm(em_df$x, mean = em_df$mu, sd = em_df$sigma)
  min <- which.min(mm$mu)

  # define threshold
  threshold <- round(mm$mu[min] + 2.5*mm$sigma[min], 3)
  ui_info("Threshold: {threshold}")
  if (threshold >= 3) {threshold <- 3}

  # define clusters
  mnp_clusters <- as.character(mnp_scores %>%
                                 group_by(harmony_leiden) %>%
                                 summarize(median_mnp_score = median(mnp_diff_score)) %>%
                                 filter(median_mnp_score > threshold & median_mnp_score > 0) %>%
                                 pull(harmony_leiden))
  nonmnp_clusters <- setdiff(unique(object$harmony_leiden), mnp_clusters)
  cc_props <- as.matrix(table(object$harmony_leiden, object$cc_diff_label))
  cc_props <- data.frame(cluster = rownames(cc_props), cycling = cc_props[, 1, drop = TRUE], noncycling = cc_props[, 2, drop = TRUE])
  cc_props <- cc_props %>% group_by(cluster) %>% mutate(props = cycling/(cycling+noncycling))
  cycling_cluster <- as.character(cc_props %>% ungroup %>% filter(props >= 0.25) %>% pull(cluster))
  nonmnp_clusters <- setdiff(nonmnp_clusters, cycling_cluster)
  mnp_clusters <- setdiff(mnp_clusters, cycling_cluster)
  object$mnp_score_label <- "Non-MNP"
  object$mnp_score_label[object$harmony_leiden %in% mnp_clusters] <- "MNP"
  object$mnp_score_label[object$harmony_leiden %in% cycling_cluster] <- "Cycling"

  # classification by voting
  counts <- object@meta.data %>%
    select(harmony_leiden, trkr20_class, trkr20_type) %>%
    group_by(harmony_leiden) %>%
    mutate(total = n()) %>%
    group_by(harmony_leiden, trkr20_type) %>%
    mutate(n = n(), prop = n/total) %>%
    distinct()
  pos_labels <- counts %>%
    filter(prop >= 0.40) %>%
    group_by(harmony_leiden) %>%
    top_n(1, prop)
  counts <- counts %>%
    filter(!harmony_leiden %in% pos_labels$harmony_leiden) %>%
    group_by(harmony_leiden) %>% top_n(2, prop) %>%
    mutate(classes = length(unique(trkr20_class))) %>%
    filter(classes == 1) %>% top_n(1, prop)
  labels <- c(pos_labels$trkr20_type, counts$trkr20_type)
  names(labels) <- c(pos_labels$harmony_leiden, counts$harmony_leiden)
  missing <- unique(object$harmony_leiden)[!unique(object$harmony_leiden) %in% names(labels)]
  if (!is_empty(missing)) {
    mixed <- rep("mixed", length(missing))
    names(mixed) <- missing
    labels <- c(labels, mixed)
  }
  labels[cycling_cluster] <- "cycling"

  # rename idents
  Idents(object) <- "harmony_leiden"
  object <- RenameIdents(object, labels)
  object[["jmp_celltype"]] <- Idents(object)

  object$mnp <- FALSE
  object$mnp[object$jmp_celltype %in% c("monocyte", "dc", "macrophage")] <- TRUE

  # label cycling cells
  if (!is_empty(cycling_cluster)) {
    cc_scores <- object[[]] %>% select(cell_barcode, harmony_leiden, mnp_diff_score) %>% filter(harmony_leiden %in% cycling_cluster)
    mnp_cycling <- cc_scores %>%
      group_by(harmony_leiden) %>%
      filter(mnp_diff_score > threshold) %>% pull(cell_barcode)
    object$mnp[mnp_cycling] <- TRUE
  }

  return(object)
}

ReprocessMNPs <- function(
  object
) {

  auc_threshold <- 0.8
  mnp_clusters <- unique(object$jmp_celltype[object$mnp == TRUE])
  nonmnp_clusters <- unique(object$jmp_celltype[object$mnp == FALSE])
  object@misc$global_markers <- presto::wilcoxauc(object, "jmp_celltype") %>% dplyr::filter(padj <= 1E-3 & auc >= 0.5 & logFC >= log(1.1))
  remove_genes <- object@misc$global_markers %>%
    group_by(group) %>%
    mutate(set = set) %>%
    filter(logFC > log(1.5) & auc > auc_threshold & padj < 1E-5) %>%
    filter(group %in% nonmnp_clusters) %>%
    arrange(desc(logFC))

  ui_done("Subsetting {sum(object$mnp)} MNPs")
  sub_object <- object[, object$mnp]
  sub_object <- DietSeurat(sub_object)
  Clean()

  sub_object <- ReduceDims(seurat.obj = sub_object, n.var.features = 3000, n.pcs = 50, remove.genes = remove_genes$feature)
  sub_object <- harmony::RunHarmony(sub_object, group.by.vars = "batch", dims.use = 1:30,
                                    max.iter.harmony = 100, verbose = TRUE, assay.use = "RNA", plot_convergence = FALSE,
                                    theta = 1, sigma = 0.1, lambda = 1, max.iter.cluster = 30)
  sub_object <- CalcPCHeuristics(sub_object, reduction = "harmony", store.name = "harmony_metrics", force_tp = TRUE)
  sub_object <- Embed(seurat.obj = sub_object, reduction = "harmony", heuristics.name = "harmony_metrics",
                      dims.use = 30, knn.use = 20)
  sub_object <- Walktrap(sub_object, "harmony_walktrap", "harmony_snn_20")
  sub_object <- AutoCluster(seurat.obj = sub_object, col.name = "harmony_leiden",
                            graph.name = "harmony_snn_20", mod.percentile = 1)
  return(sub_object)

}

VisualizeObject <- function(
  object, sub_object,
  save = TRUE
) {

  a <- DimPlot(object, reduction = "harmony_umap", group.by = "jmp_celltype",
               label = TRUE, repel = TRUE, shuffle = TRUE, label.size = 4,) +
    theme_classic(base_size = 18) +
    #ggthemes::scale_color_tableau("Tableau 20") +
    colorspace::scale_color_discrete_qualitative("Dark 3") +
    labs(x = "", y = "UMAP2", title = "") +
    RemoveAxes() + RemoveBackgrounds(outline = T) + NoLegend()
  plot(a)
  if (save) SavePlot(plot = a, filename = glue("{set}_umap_jmpcelltype"), save.data = FALSE, root = "", h = 6, w = 6)


  b <- DimPlot(object, reduction = "harmony_umap", group.by = "mnp",
               label = F, repel = F, shuffle = TRUE, label.size = 4,) +
    theme_classic(base_size = 18) +
    scale_color_manual(values = c("gray80", "gray20")) +
    labs(x = "UMAP1", y = "", title = "") +
    RemoveAxes() + RemoveBackgrounds(outline = T) + NoLegend()
  plot(b)
  if (save) SavePlot(plot = b, filename = glue("{set}_umap_mnp"), save.data = FALSE, root = "", h = 6, w = 6)

  sub_object$mnp_type_l2_plot <- gsub("_", " ", sub_object$mnp_type_l2)
  types <- levels(as.factor(sub_object$mnp_type_l2_plot))
  cols <- colorspace::qualitative_hcl(length(types), palette = "Dark 3")
  names(cols) <- types
  c <- DimPlot(sub_object, reduction = "harmony_umap", group.by = "mnp_type_l2_plot",
               label = TRUE, repel = TRUE, shuffle = TRUE, label.size = 4,) +
    theme_classic(base_size = 18) +
    #ggthemes::scale_color_tableau("Tableau 20") +
    colorspace::scale_color_discrete_qualitative("Dark 3") +
    labs(x = "", y = "", title = "") +
    RemoveAxes() + RemoveBackgrounds(outline = T) + NoLegend()
  plot(c)
  if (save) SavePlot(plot = c, filename = glue("{set}_umap_mnptypel2"), save.data = FALSE, root = "", h = 6, w = 6)

  grid <- a  + b + c + plot_layout(widths = c(1, 1, 1))
  if (save) SavePlot(plot = grid, filename = glue("{set}_umap_grid"), root = "", save.data = F, h = 4, w = 10)

  d <- ClusteredAUROCDotPlot(sub_object, "mnp_type_l2")
  if (save) SavePlot(plot = d, filename = glue("{set}_mnpmarkers"), root = "", save.data = T, h = 4, w = 6, s = 2)

  e <- ClusteredAUROCDotPlot(object, "jmp_celltype")
  if (save) SavePlot(plot = e, filename = glue("{set}_celltypemarkers"), root = "", save.data = T, h = 4, w = 6, s = 2)

}


#' Score and classify Travaglini genes
#'
#' @param object Seurat object
#' @param grouping.var Cluster variable
#'
#' @return
#' @export
TravagliniScores <- function(object, grouping.var = "seurat_clusters") {

  tra_type <- readRDS("data/reference/tra20_type.rds")
  tra_degs <- readRDS("data/reference/tra20_degs.rds")
  tra_degs <- tra_degs[tra_type$score_name]

  in_genes <- lapply(tra_degs, function(x) intersect(x, rownames(object)))
  object@meta.data[, grep("_score", colnames(object@meta.data))] <- NULL
  object <- AddModuleScore(object, features = in_genes, nbin = 20, ctrl = 30,
                           seed = 1, name = "tra20scores", search = FALSE)
  colnames(object@meta.data)[grep("tra20scores", colnames(object@meta.data))] <- names(in_genes)

  max_scores <- apply(object[[]][, names(in_genes), drop = TRUE], 1, which.max)
  max_scores <- ConvertNamedVecToDF(max_scores)
  colnames(max_scores) <- c("cell_barcode", "index")
  max_scores <- merge(max_scores, tra_type, by = "index", sort = FALSE)
  max_scores <- max_scores %>% dplyr::select(cell_barcode, score_name, class, subclass, type)
  colnames(max_scores) <- c("cell_barcode", "tra20_score_name", "tra20_class", "tra20_subclass", "tra20_type")
  rownames(max_scores) <- max_scores$cell_barcode
  max_scores <- max_scores[, -1]
  object <- AddMetaData(object, max_scores)
  head(object[[]])

  cluster_barcodes <- object[[]] %>% dplyr::select(cell_barcode, !!sym(grouping.var))
  cluster_assignments <- object[[]] %>% dplyr::select(!!sym(grouping.var), tra20_subclass)
  cluster_assignments <- cluster_assignments %>%
    group_by(!!sym(grouping.var), tra20_subclass) %>% summarize(n = n())
  cluster_assignments <- cluster_assignments %>%
    filter(rank(rev(n), ties.method = "random") == 1)
  cluster_assignments <- merge(cluster_barcodes, cluster_assignments, .by = grouping.var, sort = FALSE)
  rownames(cluster_assignments) <- cluster_assignments$cell_barcode
  cluster_assignments  <- cluster_assignments [, -1]
  object <- AddMetaData(object, cluster_assignments[, "tra20_subclass", drop = FALSE], col.name = "tra20_assignment")

  class_barcodes <- object[[]] %>% dplyr::select(cell_barcode, !!sym(grouping.var))
  class_assignments <- object[[]] %>% dplyr::select(!!sym(grouping.var), tra20_class)
  class_assignments <- class_assignments %>%
    group_by(!!sym(grouping.var), tra20_class) %>% summarize(n = n()) %>%
    filter(rank(rev(n), ties.method = "random") == 1)
  class_assignments <- merge(cluster_barcodes, class_assignments, .by = grouping.var, sort = FALSE)
  rownames(class_assignments) <- class_assignments$cell_barcode
  class_assignments  <- class_assignments [, -1]
  object <- AddMetaData(object, class_assignments[, "tra20_class", drop = FALSE], col.name = "tra20_lineage")
  return(object)

}

ClassifyMNP <- function(object) {
  tra_type <- readRDS("data/reference/tra20_type.rds")
  nonmnp <- as.character(tra_type$score_name[tra_type$mnp == 0])
  nonmnp <- nonmnp[nonmnp %in% colnames(object@meta.data)]
  mnp <- as.character(tra_type$score_name[tra_type$mnp == 1])
  mnp_diff_score <- rowSums(object@meta.data[, mnp]) - rowSums(object@meta.data[, nonmnp])
  object$mnp_diff_score <- mnp_diff_score

  mnp_scores <- object[[]] %>% select(harmony_leiden, mnp_diff_score)
  mm <- mixtools::normalmixEM(mnp_scores$mnp_diff_score, k = 2)
  x_data <- with(mm, seq(min(x), max(x), len = 1000))
  pars <- with(mm, data.frame(comp = colnames(posterior), mu, sigma, lambda))
  em_df <- data.frame(x = rep(x_data, each = nrow(pars)), pars)
  em_df$y <- em_df$lambda * dnorm(em_df$x, mean = em_df$mu, sd = em_df$sigma)
  min <- which.min(mm$mu)

  threshold <- mm$mu[min] + 2.5*mm$sigma[min]
  ui_info("Threshold: {threshold}")

  if (threshold >= 3) {threshold <- 3}
  mnp_clusters <- as.character(mnp_scores %>%
                                 group_by(harmony_leiden) %>%
                                 summarize(median_mnp_score = median(mnp_diff_score)) %>%
                                 filter(median_mnp_score > threshold & median_mnp_score > 0) %>%
                                 pull(harmony_leiden))
  nonmnp_clusters <- setdiff(unique(object$harmony_leiden), mnp_clusters)

  cc_props <- as.matrix(table(object$harmony_leiden, object$cc_diff_label))
  cc_props <- data.frame(cluster = rownames(cc_props), cycling = cc_props[, 1, drop = TRUE], noncycling = cc_props[, 2, drop = TRUE])
  cc_props <- cc_props %>% group_by(cluster) %>% mutate(props = cycling/(cycling+noncycling))
  cycling_cluster <- as.character(cc_props %>% ungroup %>% filter(props >= 0.25) %>% pull(cluster))

  nonmnp_clusters <- setdiff(nonmnp_clusters, cycling_cluster)
  mnp_clusters <- setdiff(mnp_clusters, cycling_cluster)

  counts <- object@meta.data %>%
    select(harmony_leiden, tra20_type) %>%
    filter(harmony_leiden %in% mnp_clusters) %>%
    group_by(harmony_leiden, tra20_type) %>%
    summarize(n = n()) %>%
    mutate(prop = n/sum(n))# %>%
  #filter(prop > 0.1 & prop < 0.50)
  contam_clusters <- unique(as.character(counts$harmony_leiden[!counts$tra20_type %in% c("monocyte", "dc", "macrophage")]))
  counts <- object@meta.data %>%
    select(harmony_leiden, tra20_type) %>%
    filter(harmony_leiden %in% mnp_clusters) %>%
    group_by(harmony_leiden, tra20_type) %>%
    summarize(n = n()) %>%
    mutate(prop = n/sum(n)) %>%
    filter(prop >= 0.50)
  highcontam_clusters <- unique(as.character(counts$harmony_leiden[!counts$tra20_type %in% c("monocyte", "dc", "macrophage")]))
  mnp_clusters <- setdiff(mnp_clusters, contam_clusters)
  mnp_clusters <- setdiff(mnp_clusters, highcontam_clusters)

  object$base_labels <- NA
  object$base_labels[object$harmony_leiden %in% mnp_clusters] <- "mnp"
  object$base_labels[object$harmony_leiden %in% contam_clusters] <- "prelim_mnp"
  object$base_labels[object$harmony_leiden %in% c(nonmnp_clusters, highcontam_clusters)] <- "non_mnp"
  object$base_labels[object$harmony_leiden %in% cycling_cluster] <- "cycling"

  object$mnp <- FALSE
  object$mnp[object$harmony_leiden %in% c(mnp_clusters, contam_clusters)] <- TRUE

  if (!is_empty(cycling_cluster)) {
    cc_scores <- object[[]] %>% select(cell_barcode, harmony_leiden, mnp_diff_score) %>% filter(harmony_leiden == cycling_cluster)
    mnp_cycling <- cc_scores %>%
      group_by(harmony_leiden) %>%
      filter(mnp_diff_score > threshold) %>% pull(cell_barcode)
    object$mnp[mnp_cycling] <- TRUE
  }

  ui_info("Found {sum(object$mnp)} MNP cells")
  return(object)

}

TransferRefLabels <- function(object, ref) {
  object$ref_directlabels <- ref$mps_celltype[Cells(object)]
  anchors <- FindTransferAnchors(reference = ref, query = object, reference.reduction = "pca",
                                 normalization.method = "LogNormalize", reference.assay = "integrated", n.trees = 20,
                                 k.filter = 100, k.score = 20, reduction = "pcaproject", npcs = 30, dims = 1:30)
  object <- TransferData(anchorset = anchors, reference = ref, query = object,
                         refdata = list(ref_transferlabels = "mps_celltype"))
  Idents(object) <- "mnp_type_l2"
  idents <- levels(Idents(object))
  cluster_props <- object[[]] %>%
    dplyr::group_by(mnp_type_l2, predicted.ref_transferlabels) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    mutate(mnp_type_l2 = factor(mnp_type_l2, levels = idents)) %>%
    arrange(mnp_type_l2, desc(prop)) %>%
    filter(prop >= 0.05)
  cluster_props <- cluster_props %>% unite("comb", mnp_type_l2, predicted.ref_transferlabels, remove = FALSE)

  label_props <- object[[]] %>%
    dplyr::group_by(predicted.ref_transferlabels, mnp_type_l2) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    filter(prop >= 0.05)

  # if a cluster is majority labeled, label it that predicted label
  maj_labels_df <- cluster_props %>% group_by(mnp_type_l2) %>% filter(prop >= 0.6) %>% top_n(1, prop)
  maj_labels <- maj_labels_df$predicted.ref_transferlabels
  names(maj_labels) <- maj_labels_df$mnp_type_l2
  object$ref_labels <- recode(Idents(object), !!!maj_labels)
  DimPlot(object, group.by = "ref_labels", label = TRUE, repel = TRUE) + NoLegend()

  # for remaining clusters, first limit to specific labeling
  label_df <- label_props %>% filter(prop >= 0.4 & predicted.ref_transferlabels != "LQC")
  label_df <- label_df %>% unite("comb", mnp_type_l2, predicted.ref_transferlabels, remove = FALSE)
  fil_cluster_props <- cluster_props %>% filter(comb %in% label_df$comb & !mnp_type_l2 %in% names(maj_labels)) %>% top_n(1, prop)
  labels_2 <- fil_cluster_props$predicted.ref_transferlabels
  names(labels_2) <- fil_cluster_props$mnp_type_l2

  labels <- c(maj_labels, labels_2)
  remaining <- cluster_props %>% filter(!mnp_type_l2 %in% names(labels)) %>% top_n(1)
  rls <- remaining$predicted.ref_transferlabels
  names(rls) <- remaining$mnp_type_l2
  labels <- c(labels, rls)
  labels[grepl("LQC", names(labels))] <- "LQC"

  object$ref_labels <- as.character(recode(Idents(object), !!!labels))
  sl_df <- label_props %>% filter(prop >= 0.6 &
                                    predicted.ref_transferlabels != "LQC" &
                                    !predicted.ref_transferlabels %in% labels &
                                    !grepl("LQC", mnp_type_l2))
  for (i in 1:nrow(sl_df)) {
    ui_info("Adding {sl_df[i, 'predicted.ref_transferlabels']}")
    cl <- sl_df[i, "mnp_type_l2", drop = T]
    lbl <- sl_df[i, "predicted.ref_transferlabels", drop = T]
    boolean <- (object$mnp_type_l2 == cl) & (object$predicted.ref_transferlabels == lbl)
    object$ref_labels[boolean] <- as.character(sl_df[i, "predicted.ref_transferlabels", drop = T])
  }
  return(object)
}
