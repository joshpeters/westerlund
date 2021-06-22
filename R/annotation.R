
# requires rentrez
# devtools::install_github("ropensci/rentrez")

ProbeNCBI <- function(gene, species) {

  if (species == "human"){
    species <- "Homo sapiens"
  } else if (species == "mouse") {
    species <- "Mus musculus"
  } else {
    stop("Check species input")
  }

  rentrez::entrez_dbs()
  rentrez::entrez_db_searchable(db = "gene")

  results <- rentrez::entrez_search(db = "gene", term = "(TGFB1[Gene Name]) AND Homo sapiens[Organism] ")

  stopifnot(length(results$count) > 0)

  info <- rentrez::entrez_summary(db = "gene", id = results$ids[1], rettype = "xml", parsed = TRUE)
  info
  info <- XML::xmlToDataFrame(info)
}