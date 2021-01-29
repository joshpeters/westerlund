#' Set Future plan to multisession
#'
#' @param num.workers Number of sessions
#' @param memory.per.worker Memory per session
#'
#' @return
#' @export
StartFutureLapply <- function(num.workers = 8, memory.per.worker = 1000) {
  future::plan(multisession, workers = num.workers, gc = FALSE)
  options(future.globals.maxSize = memory.per.worker * 1024^2)
}

#' Set Future plan to sequential
#'
#' @return
#' @export
StopFutureLapply <- function() {
  future::plan(sequential)
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


#' RemoveAxes
#'
#' Modified from Seurat::NoAxes()
#'
#' @return
#' @export
#'
#' @examples
RemoveAxes <- function (..., keep.text = FALSE, keep.ticks = FALSE)
{
  blank <- element_blank()
  no_axes_theme <- theme(axis.line.x = blank,
                         axis.line.y = blank,
                         axis.ticks.x = blank,
                         axis.ticks.y = blank,
                         axis.text.x = blank,
                         axis.text.y = blank,
                         validate = TRUE, ...)
  return(no_axes_theme)
}

RemoveBackgrounds <- function(outline = FALSE, ...)
{
  if (outline) {
    no_bg_theme <- theme(panel.background = element_rect(fill = "transparent", color = "black", size = 1),
                         plot.background = element_rect(fill = "transparent", color = "transparent", size = 0),
                         legend.background = element_rect(fill = "transparent", color = "transparent", size = 0),
                         legend.box.background = element_rect(fill = "transparent", color = "transparent", size = 0),
                         panel.grid = element_blank(),
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

SpaceAxesTitles <- function(...) {
  axes_space_theme <- theme(axis.title.x = element_text(margin = margin(16, 0, 0, 0)),
                            axis.title.y = element_text(margin = margin(0, 16, 0, 0)),
                            validate = TRUE, ...)
  return(axes_space_theme)
}

GeneralTheme <- function(base_size, ...) {
  theme_classic(base_size = base_size, ...) +
    ggeasy::easy_all_text_color("black") +
    theme (

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
