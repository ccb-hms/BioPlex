############################################################
# 
# author: Ludwig Geistlinger
# date: 2021-02-16 13:09:57
# 
# descr: obtain and clean transcriptome data from GEO
#        for HEK293 and HCT116 cells
# 
############################################################

#' @import SummarizedExperiment
#' @importFrom methods as
#' @export
getGSE122425 <- function(cache = TRUE)
{
    rname <- "GSE122425"
    # should a cached version be used
    if(cache)
    {
        se <- .getResourceFromCache(rname)
        if(!is.null(se)) return(se)
    }

    # pull from GEO
    eset <- GEOquery::getGEO(rname)[[1]]
    se <- as(eset, "SummarizedExperiment")
    
    # expression data not well formated, ...
    # we need to pull the supplementary files to get the counts
    # NOTE: .xls extension is errorneous, it's actually a .tsv file
    GEOquery::getGEOSuppFiles(rname)
    count.file <- "all.counts.293_vs_293NK.edgeR_all.xls.gz"
    count.file.prefix <- file.path(rname, rname)
    count.file <- paste(count.file.prefix, count.file, sep = "_")
    cont <- read.delim(count.file)  

    # we have raw counts and rpkms here in one matrix, ...
    # let's pull them out separately and make each one an assay
    ind <- rep("HEK293", 6)
    ind[4:6] <- paste0(ind[4:6], "NK")
    ind <- paste(ind, "SEQ", sep = ".")
    ind <- paste0(ind, rep(1:3, 2))    

    raw <- cont[,ind]    
    ind <- paste(ind, "RPKM", sep = "_")
    rpkm <- cont[,ind]    

    colnames(raw) <- colnames(rpkm) <- colnames(se)
    rownames(raw) <- rownames(rpkm) <- cont$gene_id

    # store in an SE
    cdat <- colData(se)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(raw = raw,
                                                     rpkm = rpkm),
                                                     colData = cdat)
    anno.cols <- c("GeneSymbol", "KO", "GO", "length")
    rowData(se) <- cont[,anno.cols] 

    # cache and clean up
    .cacheResource(se, rname)
    rm.files <- list.files(rname, full.names = TRUE, all.files = TRUE)
    file.remove(rm.files[-c(1,2)])
    file.remove(rname)

    return(se)
}
