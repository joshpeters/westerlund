#' Set Future plan to multisession
#'
#' @param num.workers Number of sessions
#' @param memory.per.worker Memory per session
#'
#' @return
#' @export
StartFutureLapply <- function(num.workers = 8, memory.per.worker = 1000) {
  future::plan("multisession", workers = num.workers, gc = FALSE)
  options(future.globals.maxSize = memory.per.worker * 1024^2)
}

#' Set Future plan to sequential
#'
#' @return
#' @export
StopFutureLapply <- function() {
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
