############################################################
# 
# author: Ludwig Geistlinger

# descr: misc utils
# 
############################################################

.cacheResource2 <- function(rname, fpath, ucdir = "BioPlex", ...)
{
    cache.dir <- tools::R_user_dir(ucdir, which = "cache")
    bfc <- BiocFileCache::BiocFileCache(cache.dir)
    
    # check if fpath is being tracked 
    qgsc <-  BiocFileCache::bfcquery(bfc, fpath)
    if(BiocFileCache::bfccount(qgsc) == 0) 
    {
        rpath <- BiocFileCache::bfcadd(bfc, rname, fpath, ...) 
    }
    else rpath <- qgsc$rpath
    return(rpath)
}

.getResourceFromCache2 <- function(rname, 
                                   update = TRUE,
                                   ucdir = "BioPlex",
                                   ...)
{
    cache.dir <- tools::R_user_dir(ucdir, which = "cache")
    bfc <- BiocFileCache::BiocFileCache(cache.dir)
    
    qgsc <- BiocFileCache::bfcquery(bfc, rname, exact = TRUE)

    # is there already a cached version?
    res <- NULL
    if(BiocFileCache::bfccount(qgsc))
    {
        rid <- qgsc$rid
        
        # is the cached version outdated?
        nu <- FALSE
        if(update) nu <- BiocFileCache::bfcneedsupdate(bfc, rid)
        if(!isFALSE(nu)) BiocFileCache::bfcdownload(bfc, rid, ask = FALSE, ...)
        message("Using cached version from ", qgsc$create_time)
        res <- BiocFileCache::bfcrpath(bfc, rname)
    }
    return(res)   
}
