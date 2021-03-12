############################################################
# 
# author: Ludwig Geistlinger
# date: 2021-02-26 08:22:09
# 
# descr: obtain BioPlex PPIs
# 
############################################################

#' @title Obtain BioPlex protein-protein interaction data
#' @description Functionality for retrieving the BioPlex protein-protein
#' interaction data. Available networks include:
#' \itemize{
#' \item BioPlex 293T cells (versions 1.0, 2.0, and 3.0), 
#' \item BioPlex HCT116 cells (version 1.0).
#' } See references.
#' @param cell.line character. Valid options include:
#' \itemize{ 
#' \item \code{"293T"}: derivative of human embryonic kidney 293 cell line,
#' \item \code{"HCT116"}: human colon cancer cell line 116.
#' }
#' Defaults to \code{"293T"}.
#' @param version character. Valid options include \code{"1.0"}, \code{"2.0"},
#' and \code{"3.0"} for 293T cells. For HCT116 cells, only \code{"1.0"} is 
#' available. 
#' Defaults to \code{"3.0"}.
#' @param cache logical. Should a locally cached version used if available?
#' Defaults to \code{TRUE}.
#' @return A \code{data.frame}. 
#' @references BioPlex: \url{https://bioplex.hms.harvard.edu/interactions.php}
#' @examples
#' # (1) Obtain the latest version of the 293T PPI network
#' bp.293t <- getBioPlex(cell.line = "293T", version = "3.0")
#' 
#' # (2) Obtain the latest version of the HCT116 PPI network
#' bp.hct116 <- getBioPlex(cell.line = "HCT116", version = "1.0")
#' @importFrom utils read.delim
#' @export
getBioPlex <- function(cell.line = c("293T", "HCT116"),
                       version = c("3.0", "1.0", "2.0"),
                       cache = TRUE)
{
    bioplex.url <- "https://bioplex.hms.harvard.edu/data/BioPlex_"
  
    cell.line <- match.arg(cell.line)
    version <- match.arg(version)
    
    # we only have version 1.0 for HCT116 cells currently
    if(cell.line == "HCT116") version <- "1.0"
    
    clver <- paste(cell.line, version, sep = ".")
    rname <- paste("bioplex", clver, sep = ".")
    
    # should a cache version be used?
    if(cache)
    {
      bioplex <- .getResourceFromCache(rname)
      if(!is.null(bioplex)) return(bioplex)        
    }
    
    # get the data
    file.ext <- switch(clver,
                       `293T.1.0` = "interactionList_v2",
                       `293T.2.0` = "interactionList_v4a",
                       `293T.3.0` = "293T_Network_10K_Dec_2019",
                       `HCT116.1.0` = "HCT116_Network_5.5K_Dec_2019")
    file.ext <- paste(file.ext, "tsv", sep = ".")
    ppi.file <- paste0(bioplex.url, file.ext)
    bioplex <- read.delim(ppi.file)
    
    # clean up & cache
    .cacheResource(bioplex, rname)
    return(bioplex)
}
