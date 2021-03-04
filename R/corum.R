############################################################
# 
# author: Ludwig Geistlinger
# date: 2021-02-16 09:20:02
# 
# descr: download CORUM data 
# 
############################################################

#' @importFrom methods as
#' @importFrom utils download.file unzip
#' @export
getCorum <- function(set = c("all", "core", "splice"),
                     cache = TRUE)
{
    corum.url <- "http://mips.helmholtz-muenchen.de/corum/download"
    set <- match.arg(set)
    rname <- paste("corum", set, sep = "-")

    # should a cache version be used?
    if(cache)
    {
        corum <- getResourceFromCache(rname)
        if(!is.null(corum)) return(corum)        
    }

    # download, unzip, and read in
    set <- paste0(set, "Complexes.txt.zip")
    download.file(file.path(corum.url, set), destfile = set)
    unzip(set)
    set <- sub(".zip$", "", set)
    corum <- vroom::vroom(set)

    # clean up & cache
    cacheResource(corum, rname)
    file.remove(c(set, paste0(set, ".zip")))
    return(corum) 
}

#' @export
corum2list <- function(corum.df, 
                       organism = "Human",
                       subunit.id.type = c("UNIPROT", "ENTREZID"))
{
    subunit.id.type <- match.arg(subunit.id.type)
    stopifnot(organism %in% corum.df$Organism)
    
    corum.df <- subset(corum.df, Organism == organism)
    subunit.id.type <- ifelse(subunit.id.type == "UNIPROT", "UniProt", "Entrez")
    id.col <- paste0("subunits(", subunit.id.type, " IDs)")
    
    gs <- strsplit(corum.df[[id.col]], " ?; ?")
    complex.name <- gsub(" ", "_", corum.df$ComplexName)
    names(gs) <- paste0("CORUM", corum.df$ComplexID, "_", complex.name)
    return(gs)
}
    
corum2GeneSetCollection <- function(corum.df)
{
    
}

