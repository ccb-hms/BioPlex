############################################################
# 
# author: Ludwig Geistlinger
# date: 2021-03-31 14:53:50
# 
# descr: utilities for obtaining & analysing proteome data 
# 
############################################################

#' @title Convenient access to the CCLE proteome data
#' @description Functionality for storing the protein expression data
#' from the Cancer Cell Line Encyclopedia (CCLE) in a 
#' \code{\linkS4class{SummarizedExperiment}}.
#' @param df a \code{data.frame} storing the CCLE protein expression data with
#' one measurement in each row. Typically obtained from \code{ExperimentHub}.
#' See examples. 
#' @param cell.line character. One or more cell line IDs such as \code{"HCT116"}
#' (human colon cancer cell line 116). Use \code{NULL} to not subset by cell line.
#' Defaults to \code{"HCT116"}, which will then subset the df to measurements for
#' HCT116 only.
#' @return A \code{\linkS4class{SummarizedExperiment}} storing protein expression
#' data for the specified cell line(s). 
#' @references CCLE proteomics: \url{https://gygi.hms.harvard.edu/publications/ccle.html}
#' @examples
#' 
#'   # Connect to ExperimentHub
#'   eh <- ExperimentHub::ExperimentHub()
#'
#'   # Obtain CCLE proteome data frame
#'   AnnotationHub::query(eh, c("gygi", "depmap"))
#'   ccle.prot <- eh[["EH3459"]]
#'   ccle.prot <- as.data.frame(ccle.prot)
#'
#'   # Turn into a SummarizedExperiment
#'   se <- ccleProteome2SummarizedExperiment(ccle.prot)
#'   
#' @export
ccleProteome2SummarizedExperiment <- function(df, cell.line = "HCT116")
{
    ids <- strsplit(df[,"cell_line"], "_")
    ids <- vapply(ids, `[`, character(1), x = 1)
  
    # subset by specific cell lines?
    if(!is.null(cell.line))
    {
        ind <- ids %in% cell.line
        sdiff <- setdiff(cell.line, unique(ids))
        if(length(sdiff)) warning("Ignoring invalid cell lines ", 
                                  paste(sdiff, collapse = ", "))
        df <- df[ind,]
        ids <- ids[ind]
    }

    # split + extract by cell line
    df.spl <- split(df, ids)
    
    # first: some sanity checks on dimensions and UNIPROT IDs
    nrows <- vapply(df.spl, nrow, integer(1))
    ind <- nrows == nrows[1]
    if(!all(ind)) df.spl <- df.spl[ind]
    unip <- vapply(df.spl, 
                   function(x) x[,"uniprot_acc"],
                   character(nrows[1]))
    uids <- apply(unip, 1, unique)
    stopifnot(is.character(uids))
    
    # second: pull out expression and metadata separately 
    pexpr <- vapply(df.spl, 
                    function(x) x[,"protein_expression"],
                    numeric(nrows[1]))
    rownames(pexpr) <- uids
    
    # store in SE
    se <- SummarizedExperiment(assays = list(expr = pexpr))
    rel.cols <- c("gene_name", "entrez_id")
    rowData(se) <- df.spl[[1]][,rel.cols]
    colnames(rowData(se)) <- c("SYMBOL", "ENTREZID")
    
    return(se)
}

#' @title Obtain BioPlex3 proteome data
#' @description Functionality for retrieving the BioPlex3 protein expression 
#' data comparing expression in the HCT116 and the 293T cell lines.
#' @param cache logical. Should a locally cached version used if available?
#' Defaults to \code{TRUE}.
#' @return A \code{\linkS4class{SummarizedExperiment}} storing protein expression
#' data for the both cell line(s) with 5 replicates each. 
#' @references BioPlex: \url{https://bioplex.hms.harvard.edu}
#' @examples
#' 
#'   se <- getBioplexProteome()
#'   
#' @importFrom utils read.csv
#' @export
getBioplexProteome <- function(cache = TRUE)
{
    bioplex.url <- file.path("https://bioplex.hms.harvard.edu/data", 
                             "293T_HCT116_ProteomeComparison.tsv")

    # should a cached version be used?
    rname <- "bp.prot"
    if(cache)
    {   
        se <- .getResourceFromCache(rname)
        if(!is.null(se)) return(se)
    }   
  
    # read and extract the data
    dat <- read.delim(bioplex.url)
    ind <- grep("scaled$", colnames(dat))
    emat <- dat[,ind]
    
    ids <- strsplit(dat[,"Protein.Id"], "\\|")
    ids <- vapply(ids, `[`, character(1), x = 2)
    rownames(emat) <- ids
    colnames(emat) <- sub(".scaled$", "", colnames(emat))
    emat <- as.matrix(emat)
    
    # turn into a SummarizedExperiment
    se <- SummarizedExperiment(assays = list(exprs = emat))
    se$cell.line <- rep(c("HCT116", "293T"), each = 5)
    rcols <- c("GeneID", "Gene.Symbol", "X..Peps", "Log2ratio", "AdjPValue")
    rowData(se) <- dat[,rcols]
    colnames(rowData(se)) <- c("ENTREZID", "SYMBOL", 
                               "nr.peptides", "log2ratio", "adj.pvalue")
    rowData(se)$ENTREZID <- as.character(rowData(se)$ENTREZID)
    .cacheResource(se, rname)
    return(se)
}
