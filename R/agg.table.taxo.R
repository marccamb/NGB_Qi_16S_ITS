#' Table aggregation according to taxonomy
#' Aggregates an ASV/OTU table according to a given taxonomic level
#' @param tab An abundance table containing samples in columns and OTUs/ASV in rows.
#' @param tax.table A table containing the taxonomy of each ASV/OTU.
#' @param tax.lvl A string specifying the taxonomic level to be used for aggregation.
#' \code{tax.lvl} should be the name of one colum of tax.table (case insensitive).
#' The default is "genus".
#' @return An abundance table with the sum of counts for each ASV/OTU corresponding to
#' a given taxonomic level, found in each sample. ASV/OTUs that are not assigned at this
#' specific level are discarded.
#' @import stats
#'
# 2019-11-28
# Marine Cambon

agg.table.taxo <- function(tab, tax.lvl="genus", tax.table) {
  if (any(rownames(tab) != rownames(tax.table))) stop("The ASV/OTU table and the taxonomy table should have the same
                                                      row names. Check that they are ordered the same way")
  message(paste('Table aggregation to the', tax.lvl, "level."))
  #message('Please be sure that the ASV/OTU table and the taxonomy table are ordered the same way')
  if(nrow(tab) != nrow(tax.table)) stop("The ASV/OTU table and the taxonomy table do not have the same number of rows")
  tax <- tax.table[,grep(tax.lvl, colnames(tax.table), ignore.case = T)]
  #tax[is.na(tax)] <- "Unknown"
  tab <- stats::aggregate(tab, by=list("taxo"=tax), FUN=sum)
  rownames(tab) <- tab[,1]
  tab <- tab[,-1]
  return(tab)
}
