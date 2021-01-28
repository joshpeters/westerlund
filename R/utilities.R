#' Set Future plan to multisession
#'
#' @param num.workers Number of sessions
#' @param memory.per.worker Memory per session
#'
#' @return
#' @export
StartFutureLapply <- function(num.workers = 8, memory.per.worker = 1000) {
  future::plan(strategy = "multisession", workers = num.workers, gc = FALSE)
  options(future.globals.maxSize = memory.per.worker * 1024^2)
}

#' Set Future plan to sequential
#'
#' @return
#' @export
StopFutureLapply <- function() {
  future::plan(strategy = "eager")
  usethis::ui_done("Multisession stopped.")
  gc()
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
  if(exists("geosketch")) {
    geosketch <- reticulate::import("geosketch")
  }
  stopifnot(reduction %in% names(object@reductions))
  stopifnot(ncol(object@reductions[[reduction]]) > dims)

  embeddings <- object@reductions[[reduction]]@cell.embeddings[, 1:dims]
  index <- unlist(geosketch$gs(embeddings, as.integer(num.cells)))
  sketch <- object[, index]
  return(sketch)
}
