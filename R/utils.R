############################################################
# 
# author: Ludwig Geistlinger
# date: 2021-01-05 18:08:41
# 
# descr: misc utils
# 
############################################################

.cacheResource <- function(res, rname, ucdir="BioPlex")
{
    # are we running interactively?
    cache.dir <- tools::R_user_dir(ucdir, which = "cache")
    bfc <- BiocFileCache::BiocFileCache(cache.dir)
    
    # replace existing version if necessary 
    qgsc <-  BiocFileCache::bfcquery(bfc, rname, exact = TRUE)
    if(BiocFileCache::bfccount(qgsc)) BiocFileCache::bfcremove(bfc, qgsc$rid) 
    
    cache.file <- BiocFileCache::bfcnew(bfc, rname)
    saveRDS(res, file=cache.file)
}

.getResourceFromCache <- function(rname, 
    update.value=36, update.unit="weeks", ucdir="BioPlex")
{
    # are we running interactively?
    cache.dir <- tools::R_user_dir(ucdir, which = "cache")

    bfc <- BiocFileCache::BiocFileCache(cache.dir)
    qgsc <-  BiocFileCache::bfcquery(bfc, rname, exact = TRUE)

    # is there already a cached version?
    res <- NULL
    if(BiocFileCache::bfccount(qgsc))
    {
        # is the cached version outdated?
        dt <- difftime(Sys.time(), qgsc$create_time, units=update.unit)   
        if(is.na(update.value) || dt < update.value)
        {
            message(paste("Using cached version from", qgsc$create_time))
            res <- readRDS(BiocFileCache::bfcrpath(bfc, rname))
        }
    }
    return(res)   
}
