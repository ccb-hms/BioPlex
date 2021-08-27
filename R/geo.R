############################################################
# 
# author: Ludwig Geistlinger
# date: 2021-02-16 13:09:57
# 
# descr: obtain and clean transcriptome data from GEO
#        for HEK293 and HCT116 cells
# 
############################################################

#' @title Convenient access to 293T transcriptome data from GEO
#' @description Functionality for storing the 293T RNA-seq data
#' from GSE122425 in a 
#' \code{\linkS4class{SummarizedExperiment}}. The dataset includes three
#' wild type samples and three NSUN2 knockout samples.
#' @param cache logical. Should a locally cached version used if available?
#' Defaults to \code{TRUE}.
#' @return A \code{\linkS4class{SummarizedExperiment}} storing RNA-seq
#' data for the 293T cell line. 
#' @references GSE122425: \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122425}
#' @examples
#' 
#'   # Obtain the data as a SummarizedExperiment
#'   se <- getGSE122425()
#'  
#' @import SummarizedExperiment
#' @importFrom methods as
#' @export
getGSE122425 <- function(cache = TRUE)
{
    rname <- "GSE122425"
    stub <- gsub("\\d{1,3}$", "nnn", rname, perl = TRUE)
    url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/", stub, rname)
    count.file <- paste0(rname, "_all.counts.293_vs_293NK.edgeR_all.xls.gz")
    url <- paste0(url, count.file)   

    # should a cached version be used
    if(cache) se.file <- .getResourceFromCache2(rname,
                                                update = FALSE,
                                                FUN = .gse2se)
    if(!cache || is.null(se.file))
    {
        se.file <- .cacheResource2(rname, url, download = FALSE, ext = ".rds")
        se.file <- suppressMessages(.getResourceFromCache2(rname, FUN = .gse2se))
    }
    se <- readRDS(se.file)
    return(se)
}

.gse2se <- function(from, to)
{
    # pull from GEO
    eset <- GEOquery::getGEO("GSE122425")[[1]]
    se <- as(eset, "SummarizedExperiment")
    
    # expression data not well formated, ...
    # we need to pull the supplementary files to get the counts
    # NOTE: .xls extension is errorneous, it's actually a .tsv file
    cont <- read.delim(from)  

    # we have raw counts and rpkms here in one matrix, ...
    # let's pull them out separately and make each one an assay
    ind <- rep("HEK293", 6)
    ind[4:6] <- paste0(ind[4:6], "NK")
    ind <- paste(ind, "SEQ", sep = ".")
    ind <- paste0(ind, rep(seq_len(3), 2)) 

    raw <- as.matrix(cont[,ind])    
    ind <- paste(ind, "RPKM", sep = "_")
    rpkm <- as.matrix(cont[,ind])    

    colnames(raw) <- colnames(rpkm) <- colnames(se)
    rownames(raw) <- rownames(rpkm) <- cont$gene_id

    # store in an SE
    cdat <- colData(se)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(raw = raw,
                                                     rpkm = rpkm),
                                                     colData = cdat)
    anno.cols <- c("GeneSymbol", "KO", "GO", "length")
    rowData(se) <- cont[,anno.cols] 
    colnames(rowData(se))[1] <- "SYMBOL"

    saveRDS(se, file = to)
    return(TRUE)
}

