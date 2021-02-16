############################################################
# 
# author: Ludwig Geistlinger
# date: 2021-02-16 09:20:02
# 
# descr: download CORUM data 
# 
############################################################

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
    cacheResource(corum, rname)

    # clean up
    file.remove(c(set, paste0(set, ".zip")))
    return(corum) 
}
