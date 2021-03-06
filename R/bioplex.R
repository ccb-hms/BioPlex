############################################################
# 
# author: Ludwig Geistlinger
# date: 2021-02-26 08:22:09
# 
# descr: obtain BioPlex PPIs
# 
############################################################

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
