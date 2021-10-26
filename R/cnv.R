############################################################
# 
# author: Ludwig Geistlinger
# date: 2021-02-26 08:22:09
# 
# descr: obtain CNV data for HEK293T and HCT116 cells
# 
############################################################

#' @title Obtain HEK293 genome data 
#'
#' @description Functionality for retrieving genomic data for different
#' lineages of the human embryonic kidney HEK293 cell line. 
#' Returned genomic coordinates are based on the *hg18* human genome assembly.
#' See references.
#' @param track character. Genome track to retrieve. Valid options include: 
#' \itemize{ 
#' \item \code{"cnvhmm"}: regions of copy number variation (CNV) as inferred by
#' a hidden Markov model (HMM) algorithm,
#' \item \code{"cnvsnp"}: CNV regions as inferred from Illumina SNP arrays
#' }
#' Defaults to \code{"cnvhmm"}.
#' @param cell.line character. Valid options include:
#' \itemize{ 
#' \item \code{"293T"}: highly-transfective derivative of human embryonic kidney
#' 293 cell line,
#' }
#' Defaults to \code{"293T"}.
#' @param cache logical. Should a locally cached version used if available?
#' Defaults to \code{TRUE}.
#' @return A \code{GRanges} object storing genomic coordinates and genomic scores
#' of regions of interest. 
#' @references \url{http://hek293genome.org}
#' @examples
#'    cnv.hmm <- getHEK293GenomeTrack(track = "cnv.hmm", cell.line = "293T")
#' 
#' @importFrom utils read.delim
#' @export
getHEK293GenomeTrack <- function(track = c("cnv.hmm", "cnv.snp"),
                                 cell.line = "293T",
                                 cache = TRUE)
{
    hek.url <- "http://bioinformatics.psb.ugent.be/downloads/genomeview/hek293"
    
    cell.line <- match.arg(cell.line)
    track <- match.arg(track)
    rname <- paste(cell.line, track, sep = ".")
 
    # should a cache version be used?
    if(cache) track.file <- .getResourceFromCache2(rname)
    if(!cache || is.null(track.file))    
    {   
        track <- switch(track,
            cnv.hmm = "CNV/cnv-293T-ALL.bed",
            cnv.snp = "SNP_array/T_chr_hg18.bed")

        track.file <- file.path(hek.url, track)
        track.file <- .cacheResource2(rname, track.file)
    }

    track <- read.delim(track.file, skip = 1, header = FALSE) 
    track <- track[,1:4]
    colnames(track) <- c("chr", "start", "end", "score")
    track <- GenomicRanges::makeGRangesFromDataFrame(track, keep.extra.columns = TRUE)
    GenomeInfoDb::genome(track) <- "hg18"
    return(track)
}


