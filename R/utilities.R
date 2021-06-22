
#' Convert count matrix to sparse matrix
#'
#' @param matrix
#' @param gene_column Logical if first column is genes
#'
#' @return dgCMatrix
ConvertCSVMatrix <- function(matrix, gene_column = TRUE) {
  genes <- matrix[, 1, drop = TRUE]
  barcodes <- colnames(matrix)[2:ncol(matrix)]
  matrix <- matrix[, -1]
  matrix <- Matrix::Matrix(as.matrix(matrix), sparse = TRUE)
  colnames(matrix) <- barcodes
  rownames(matrix) <- genes
  return(matrix)
}

#' Set Future plan to multisession
#'
#' @param num.workers Number of sessions
#' @param memory.per.worker Memory per session
#'
#' @return
#' @export
StartFuture <- function(num.workers = 4, memory.per.worker = 1000) {
  future::plan("multisession", workers = num.workers, gc = FALSE)
  options(future.globals.maxSize = memory.per.worker * 1024^2)
}

#' Set Future plan to sequential
#'
#' @return
#' @export
StopFuture <- function() {
  future::plan("sequential")
  usethis::ui_done("Multisession stopped.")
  Clean()
}

#' Geosketch Seurat object
#'
#' See \href{https://github.com/brianhie/geosketch}{geosketch}\cr\cr
#' Requires geosketch installation and proper reticulate setup
#'
#' @param object Seurat object
#' @param reduction Reduced dimensions to use for sketching
#' @param dims Number of dimensions to use
#' @param num.cells Number of desired cells
#'
#' @return sketch Seurat object downsampled to desired number of cells
#' @export
Geosketch <- function(object, reduction, dims, num.cells) {
  if(!exists("geosketch")) {
    geosketch <- reticulate::import("geosketch")
  }
  stopifnot(reduction %in% names(object@reductions))
  stopifnot(ncol(object@reductions[[reduction]]) >= dims)

  embeddings <- object@reductions[[reduction]]@cell.embeddings[, 1:dims]
  index <- unlist(geosketch$gs(embeddings, as.integer(num.cells)))
  sketch <- object[, index]
  return(sketch)
}

#' Timestamp
#' Get current timestamp
#' @return timestamp
#' @export
Timestamp <- function() {
  return(format(Sys.time(), '%y%m%d%I%M'))
}

#' Timestamp
#' Get current timestamp
#' @return timestamp
#' @usage GetDate()
#' @export
GetDate <- function() {
  return(format(Sys.Date(), "%y%m%d"))
}

#' Clean
#' gc but invisible
#' @return NULL
#' @export
Clean <- function() {
  return(invisible(gc(verbose = FALSE)))
}

#' Convert variable to string
#'
#' This function turns variable names into a string.
#' @param variable variable to convert to a string name
#' @return Returns a string of variable name
#' @usage VarToString(variable)
#' @export
VarToString <- function(variable) {
  return(deparse(substitute(variable)))
}

# See https://github.com/satijalab/seurat/blob/master/R/utilities.R for following functions
#' Set a default value if an object is null
#'
#' @param lhs An object to set if it's null
#' @param rhs The value to provide if x is null
#'
#' @return rhs if lhs is null, else lhs
#'
#' @author Hadley Wickham
#' @references https://adv-r.hadley.nz/functions.html#missing-arguments
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

#' Set a default value if an object is NOT null
#'
#' @param lhs An object to set if it's NOT null
#' @param rhs The value to provide if x is NOT null
#'
#' @return lhs if lhs is null, else rhs
#'
#' @author Hadley Wickham
#' @references https://adv-r.hadley.nz/functions.html#missing-arguments
`%iff%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(rhs)
  } else {
    return(lhs)
  }
}

#' Convert gene names based on Ensembl v86 and HGNC databases
#' @param genes List of genes to convert
#' @param verbose Verbosity toggle
#' @param hgnc HGNC database to use
#' @return Gene ids and conversions data frame
ConvertGeneNames <- function(genes, verbose = TRUE) {

  hgnc <- readRDS("data/reference/hgnc.rds")

  # map IDs using the most recent ENSEMBL version
  num_genes <- length(genes)
  assertthat::assert_that(length(unique(genes)) == num_genes)
  ids <- data.frame(orig_name = genes)
  ensembl_mapped <- ensembldb::mapIds(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, keys = as.character(genes), keytype = "SYMBOL", column = "GENEID")
  ensembl_mapped <- ConvertNamedVecToDF(ensembl_mapped)
  ids <- merge(ids, ensembl_mapped, by.x = "orig_name", by.y = "name", all.x = TRUE, sort = FALSE)
  colnames(ids) <- c("orig_name", "id")
  if (verbose) usethis::ui_info("{sum(duplicated(ids$id))} duplicated Ensembl mappings")
  if (verbose) usethis::ui_oops("{sum(is.na(ids$id))} ({round(sum(is.na(ids$id))/num_genes*100, 2)}%) unidentified genes after Ensembl mapping")
  assertthat::assert_that(all.equal(as.character(ids$orig_name), as.character(genes)))
  ids <- ids %>% mutate_all(as.character)

  # check HGNC alias database
  unids <- ids[is.na(ids$id), ]
  unids$lower <- tolower(unids$orig_name)
  hgnc$lower <- tolower(hgnc$symbol)
  unids <- merge(unids, hgnc, by = "lower", all.x = TRUE, sort = FALSE)
  colnames(unids) <- c("lower", "orig_name", "ensembl_id", "hgnc_symbol", "hgnc_id")
  unids <- unids[!duplicated(unids$orig_name), ]
  unids <- unids %>% mutate_all(as.character)
  unids <- unids[, c("orig_name", "hgnc_symbol", "hgnc_id")]
  ids <- merge(ids, unids, by = "orig_name", all.x = TRUE, sort = FALSE)
  assertthat::assert_that(length(unique(ids$orig_name)) == num_genes)
  ids$combined_id <- as.character(ids$id)
  ids$combined_id[is.na(ids$combined_id)] <- as.character(ids$hgnc_id[is.na(ids$combined_id)])
  if (verbose) usethis::ui_oops("{sum(is.na(ids$combined_id))} ({round(sum(is.na(ids$combined_id))/num_genes*100, 2)}%) unidentified genes remain after alias lookup")
  colnames(ids) <- c("orig_name", "ensembl_id", "hgnc_symbol", "hgnc_id", "id")

  # reidentify gene symbols
  symbols <- ensembldb::mapIds(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, keys = as.character(ids$id[!is.na(ids$id)]), keytype = "GENEID", column = "SYMBOL")
  symbols <- ConvertNamedVecToDF(symbols)
  ids <- merge(ids, symbols, by.x = "id", by.y = "name", all.x = TRUE, sort = FALSE)
  colnames(ids)[6] <- "ensembl_name"
  if (verbose) usethis::ui_done("{round(sum(ids$orig_name %in% ids$ensembl_name)/nrow(ids)*100, 2)}% agreement")
  ids <- ids[!duplicated(ids), ]
  ids <- ids %>% mutate_all(as.character)
  ids$name <- ids$ensembl_name
  ids$name[is.na(ids$name)] <- as.character(ids$orig_name[is.na(ids$name)])
  ids$name[ids$name %in% ids$name[duplicated(ids$name)] & ids$orig_name != ids$name] <- ids$orig_name[ids$name %in% ids$name[duplicated(ids$name)] & ids$orig_name != ids$name]

  assertthat::assert_that(nrow(ids) == num_genes)
  assertthat::assert_that(any(duplicated(ids$name)) == FALSE)

  rownames(ids) <- ids$orig_name
  return(ids)
}

#' Convert named vector to dataframe
#'
#' @param x Named vector
#'
#' @return data.frame
#' @export
ConvertNamedVecToDF <- function(x) {
  df <- data.frame(name = names(x), value = x)
  return(df)
}

SelfName <- function(list) {
  names(list) <- list
  return(list)
}
