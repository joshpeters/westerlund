# See Smillie, C.S., Biton, M.B., Ordovas, J.O., et al. (Cell, 2019) for more details
#
# Fits a contamination model using the specified "groups" (coarse-level clusters, e.g. Epithelial/Myeloid/Lymphoid/etc)
# Then applies this model to each cell subset in "idents" (fine-level clusters, e.g. Enterocytes/Macrophages/Tregs/etc)

#' Fit models to cell groups
#'
#' @param data TP10K matrix genes by cells
#' @param lineages cell lineage vector
#' @param types cell type vector
#' @param batches batch vector
#' @param anno named list mapping lineage to types
#' @param fit.n number of genes to dit
#'
#' @return list of models
#' @export
DetectAmbientRNA <- function(data,
                             lineages,
                             types,
                             batches,
                             anno = NULL,
                             fit.n = 50)
{

  # fit model
  usethis::ui_todo("Fitting models to lineages to define consensus model")
  lineage_models <- FitAmbientModel(data = data,
                                    idents = lineages,
                                    batches = batches,
                                    anno = NULL,
                                    fit.n = fit.n)

  # average model across groups
  coefs <- sapply(lineage_models, function(a) a$coefs)
  global_coefs <- apply(coefs, 1, median)

  # run average model on classes
  usethis::ui_done("Fitting consensus model to cell lineages")
  lineage_models <- FitAmbientModel(data = data,
                                    idents = lineages,
                                    batches = batches,
                                    anno = NULL,
                                    fit.n = fit.n,
                                    coefs = global_coefs)

  # run average model on types
  usethis::ui_done("Fitting consensus model to cell types")
  type_models <- FitAmbientModel(data = data,
                                 idents = types,
                                 batches = batches,
                                 anno = NULL,
                                 fit.n = fit.n,
                                 coefs = global_coefs)

  # return data
  return(list(lineage_models = lineage_models, type_models = type_models))

}



#' Fit model
#'
#' @param data TP10K matrix genes by cells
#' @param idents identities to fit
#' @param batches batches
#' @param anno
#' @param coefs.use model coefficients to use
#' @param fit.n number of genes to fit
#' @param fit.cutoff residual cutoff
#'
#' @return list describing model outputs
#' @export
FitAmbientModel <- function(data,
                            idents,
                            batches,
                            anno = NULL,
                            coefs.use = NULL,
                            fit.n = 50,
                            fit.cutoff = 2)
{

  # initialize variables
  if(is.null(anno)){
    anno <- structure(unique(as.character(idents)), names = unique(as.character(idents)))
  }
  if(any(! idents %in% anno)){
    stop("! idents %in% anno")
  }
  groups <- names(anno)

  # iterate over groups
  res <- pbapply::pbsapply(groups, function(group) {

    # subset data
    i <- idents %in% anno[[group]]
    j <- idents %in% setdiff(idents, i)

    # sample frequencies
    f <- table(as.factor(batches)[i])
    f <- as.matrix(f/sum(f))

    # group mean
    u1 <- Matrix::rowMeans(data[,i])

    # other means
    u2 <- sapply(unique(batches), function(a){
      Matrix::rowSums(data[,j & (batches  ==   a)])/sum(batches  ==   a)
    })

    stopifnot(colnames(u2)  ==   colnames(f))
    u2 <- (u2 %*% f)[,1]

    l1 <- Log2RmInf(u1)
    l2 <- Log2RmInf(u2)

    # fit boundaries
    lo <- quantile(l2[u1  ==   0], .9, na.rm = T)
    hi <- sort(l2, decreasing = T)[100]
    exclude <- list(c(-Inf, lo), c(hi, Inf))

    # select points for regression
    genes_fit <- names(SelectPoints(l2, l1, n  =  fit.n, dir  =  "down", nbins  =  10, loess  =  TRUE, exclude  =  exclude))
    fit <- MASS::rlm(l1[genes_fit] ~ l2[genes_fit])
    coefs <- as.matrix(coef(fit))

    if(!is.null(coefs.use)){
      coefs <- cbind(coefs.use, coefs)
    }

    # calculate residuals
    residuals <- l1 - (coefs[2,1]*l2 + coefs[1,1])
    genes_ambient <- names(which(residuals < 2))

    # update results
    list(u1  =  u1, u2  =  u2, fit  =  fit, coefs  =  coefs, residuals  =  residuals, lab.fit  =  genes_fit, lab.con  =  genes_ambient)

  }, simplify  =  F)

  names(res) <- groups
  return(res)
}

#' Plot Ambient RNA model fits
#'
#' @param u1 mean expression in-group
#' @param u2 mean expreession out-group
#' @param coefs model coefficients
#' @param residuals model residuals
#' @param genes_label genes to label
#' @param genes_fit genes used for fitting
#' @param label_n number of genes to label
#' @param residual_cutoff residual cutoff
#'
#' @return ggplot2
#' @export
PlotAmbientRNA <- function(u1,
                           u2,
                           coefs,
                           residuals,
                           genes_label = NULL,
                           genes_fit = NULL,
                           label_n = 10,
                           residual_cutoff = 2)
{

  # log-transform
  l1 <- log2(u1 + .5*min(u1[u1 > 0]))
  l2 <- log2(u2 + .5*min(u2[u2 > 0]))

  # contamination
  ambient_genes <- names(which(residuals < residual_cutoff))

  # scatterplot data
  d <- data.frame(x = l2, y = l1,
                  lab = ifelse(names(l1) %in% genes_fit, names(l1), ""),
                  Type = rep("Other", length(l1)),
                  stringsAsFactors = F)
  d[ambient_genes, "Type"] <- "Contamination"
  d[genes_fit, "Type"] <- "Fit"
  genes_label <- intersect(rownames(d), genes_label)
  d[genes_label, "Type"] <- "Label"
  d[genes_label, "lab"] <- genes_label

  # select subset to label
  if(!is.null(label_n)){
    lab.sub <- sample(genes_fit, min(label_n, length(genes_fit)))
    d[(!d$lab %in% genes_label & !d$lab %in% lab.sub), "lab"] <- ""
  }

  # make plot
  p <- ggplot(d, aes(x = x, y = y)) +
    geom_point(aes(color = Type)) +
    ggrepel::geom_text_repel(aes(label = lab), size = 4, segment.color = "grey",
                             max.overlaps = Inf, box.padding = unit(1, "lines")) +
    xlab("log2(TP10K) (non-group)") +
    ylab("log2(TP10K) (group)") +
    scale_colour_manual(values = c("darkred", "darkorange", "gray75", "gray75")) +
    theme_classic(base_size = 18) + SpaceAxesTitles() + RemoveBackgrounds(outline = TRUE)

  # add regression lines
  coefs <- as.matrix(coefs)
  for(j in 1:ncol(coefs)){
    x0 <- (min(l1, na.rm = T) - coefs[1,j] - residual_cutoff)/coefs[2,j]
    x1 <- max(l2, na.rm = T)
    d.line <- data.frame(x = c(x0, x1))
    d.line$y <- coefs[2,j]*d.line$x + coefs[1,j] + residual_cutoff
    if(j  ==   1){lty <- "longdash"} else {lty <- "dotted"}
    p <- p + geom_line(data = d.line, aes(x = x, y = y), lty = lty)
  }

  return(p)

}

#' Select points
#'
#' @param x x
#' @param y y
#' @param n n
#' @param dir direcetionality
#' @param loess use Loess
#' @param nbins number of bins
#' @param bin.type type of bin
#' @param exclude exclud
#'
#' @return points vector
#' @export
SelectPoints <- function(x, y, n,
                         dir = "both",
                         loess = FALSE,
                         nbins =  25,
                         bin.type  = "equal_width",
                         exclude =  c(0,0))
{

  # fix inputs
  if(!is.list(exclude)){exclude  =  list(exclude)}
  i <- ((is.na(x) | is.na(y)) | (is.infinite(x) | is.infinite(y)))
  i <- which(!i)
  xi <- x[i]
  yi <- y[i]

  # de-trend
  if(loess  == TRUE){
    l <- LoessRegression(yi ~ xi, family = "symmetric")
    yi <- l$fit$residuals
  }

  # exclude data
  j <- apply(sapply(exclude, function(e) (e[[1]] < xi) & (xi < e[[2]])), 1, any)
  j <- which(!j)
  xi <- xi[j]
  yi <- yi[j]
  i <- i[j]

  # bin x-axis
  if(bin.type  == "equal_width"){
    groups <- Hmisc::cut2(xi, cuts = seq(from = min(xi, na.rm = T), to = max(xi, na.rm = T), length.out = nbins), m = 2*n/nbins)
  } else {
    groups <- Hmisc::cut2(xi, g = nbins)
  }

  # points
  j <- c()

  # up points
  if(dir %in% c("up", "both")){
    j <- c(j, as.numeric(Downsample(cells = 1:length(xi), groups = groups, quality.vector = yi, total.cells = n)))
  }

  # down points
  if(dir %in% c("down", "both")){
    j <- c(j, as.numeric(Downsample(cells = 1:length(xi), groups = groups, quality.vector = -1*yi, total.cells = n)))
  }

  return(i[j])
}

#' Downsample cells evenly across groups
#'
#' @param cells vector of cell barcodes
#' @param groups vector of cell groups
#' @param quality.vector vector of cell quality
#' @param total.cells total number of cells to downsample to
#' @param cells.per.group total number of cells to downsample from each group
#'
#' @return Vector of downsampled cells
#' @export
Downsample <- function(cells,
                       groups,
                       quality.vector = NULL,
                       total.cells = NULL,
                       cells.per.group = NULL)
{


  # Set ngene (default <- 1)
  if(is.null(quality.vector)){
    quality.vector <- structure(rep(1, length(cells)), names = cells)
  }
  if(is.null(names(quality.vector))){names(quality.vector) <- cells}

  # Calculate group sizes
  groups <- as.factor(groups)
  num_cells_per_group <- NumCellsPerGroup(groups = groups, total.cells = total.cells, cells.per.group = cells.per.group)

  # Downsample cells within each group
  downsampled_cells <- sapply(levels(groups), function(a){

    # Shuffle cells within group
    cells <- Resample(cells[groups  ==   a])

    # Select by highest quality
    cells[order(quality.vector[cells], decreasing = T)[1:num_cells_per_group[[a]]]]

  })
  downsampled_cells <- as.character(na.omit(unname(unlist(downsampled_cells))))
  return(downsampled_cells)
}

#' Calculate number of cells to downsample per group
#'
#' @param groups vector of groups for each cell barcode
#' @param total.cells total number of cells
#' @param cells.per.group number of cells to select per group
#'
#' @return Number of cells per group
#' @export
NumCellsPerGroup <- function(groups,
                             total.cells = NULL,
                             cells.per.group = NULL)
{

  num_cells <- sort(table(groups))

  if(!is.null(cells.per.group)){
    num_cells[num_cells > cells.per.group] <- cells.per.group
  } else {
    n <- sort(table(groups))
    if(length(n)  ==   1){
      num_cells <- total.cells
      names(num_cells) <- names(n)
    } else {
      u <- c(0, cumsum(n)[1:(length(n)-1)])
      i <- (total.cells - u)/seq(length(n), 1, -1) < n
      if(sum(i) > 0){
        num_cells[i] <- as.integer(ceiling((total.cells - sum(n[!i]))/sum(i)))
      }
    }
  }

  return(num_cells)

}

#' Sample if length(x) != 1
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
Resample <- function(x, ...){
  if(length(x) == 1) x else sample(x, ...)
}


#' Take log2 and remove infinite results
#'
#' @param x
#'
#' @return Log2 transformed input vector, `Inf = NA`
#' @export
Log2RmInf <- function(x){
  y <- log2(x)
  y[is.infinite(y)] <- NA
  return(y)
}

#' Loess regression
#'
#' @param ...
#'Î
#' @return List of model `fit`, `x` values, and fitted `y` values
#' @export
LoessRegression <- function(...){
  fit <- loess(...)
  p.x <- fit$x[order(fit$x)]
  p.y <- fit$fitted[order(fit$x)]
  return(list(fit  =  fit, x  =  p.x, y  =  p.y))
}
