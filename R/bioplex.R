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
#' @param remap.uniprot.ids logical. Should the protein-to-gene mappings from BioPlex
#' (i.e. UNIPROT-to-SYMBOL and UNIPROT-to-ENTREZID) be updated using Bioc annotation
#' functionality?
#' Defaults to \code{FALSE} which will then keep the mappings provided by BioPlex.
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
                       remap.uniprot.ids = FALSE,
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
    if(cache) ppi.file <- .getResourceFromCache2(rname)
    if(!cache || is.null(ppi.file))        
    {
        # get the data
        file.ext <- switch(clver,
                           `293T.1.0` = "interactionList_v2",
                           `293T.2.0` = "interactionList_v4a",
                           `293T.3.0` = "293T_Network_10K_Dec_2019",
                           `HCT116.1.0` = "HCT116_Network_5.5K_Dec_2019")
        file.ext <- paste(file.ext, "tsv", sep = ".")
        ppi.file <- paste0(bioplex.url, file.ext)
        ppi.file <- .cacheResource2(rname, ppi.file)
    }
    
    bioplex <- read.delim(ppi.file)
    colnames(bioplex) <- c("GeneA", "GeneB",
                           "UniprotA", "UniprotB",
                           "SymbolA", "SymbolB",
                           "pW", "pNI", "pInt")
    bioplex$GeneA <- as.character(bioplex$GeneA)
    bioplex$GeneB <- as.character(bioplex$GeneA)
    
    # remap gene ids
    if(remap.uniprot.ids) bioplex <- .remapUniprotIdsBP(bioplex)
    
    # clean up & cache
    return(bioplex)
}

#' @title Representation of BioPlex PPIs in a graph data structure
#' @description Representation of BioPlex PPIs in a \code{graphNEL} object
#' from the \code{graph} package.
#' @param bioplex.df a \code{data.frame} storing the Bioplex PPIs in a flat
#' from-to format. Typically obtained via \code{\link{getBioPlex}}.
#' @return An object of class \code{graphNEL}. 
#' @references BioPlex: \url{https://bioplex.hms.harvard.edu/interactions.php}
#' @seealso \code{\link{getBioPlex}}, \code{\link{ftM2graphNEL}}
#' @examples
#' # (1) Obtain the latest version of the 293T PPI network
#' bp.293t <- getBioPlex(cell.line = "293T", version = "3.0")
#' 
#' # (2) Turn the data into a graph 
#' bp.gr <- bioplex2graph(bp.293t)
#' 
#' @export
bioplex2graph <- function(bioplex.df) 
{
    stopifnot(is.data.frame(bioplex.df))
    node.cols <- paste0("Uniprot", c("A", "B"))
    ftm <- as.matrix(bioplex.df[,node.cols])
    ftm <- sub("-[0-9]+$", "", ftm)
    ind <- !duplicated(ftm)
    ftm <- ftm[ind,]
    bioplex.df <- bioplex.df[ind,]
    gr <- graph::ftM2graphNEL(ftm, edgemode = "directed")
    gr <- .annotateBioplexGraph(gr, bioplex.df)
    return(gr)
}



#' @title Annotate PFAM domains to BioPlex PPI graph
#' @description This function adds PFAM domain annotations to the node metadata
#' of the BioPlex PPI graph.
#' @param bp.gr an object of class \code{\linkS4class{graph}} storing the
#' BioPlex PPIs. Typically obtained via \code{\link{bioplex2graph}}.
#' @param orgdb an \code{orgdb} object storing annotation data for human.
#' @return An object of class \code{graphNEL} containing PFAM domain annotations
#' in the \code{nodeData}. 
#' @references 
#'  BioPlex: \url{https://bioplex.hms.harvard.edu/interactions.php}
#'
#'  PFAM: \url{http://pfam.xfam.org}
#' @seealso \code{\link{nodeData}}
#' @examples
#' # (1) Obtain the latest version of the 293T PPI network
#' bp.293t <- getBioPlex(cell.line = "293T", version = "3.0")
#' 
#' # (2) Turn the data into a graph 
#' bp.gr <- bioplex2graph(bp.293t)
#' 
#' # (3) Obtain orgdb package from AnnotationHub
#' ah <- AnnotationHub::AnnotationHub()
#' orgdb <- AnnotationHub::query(ah, c("orgDb", "Homo sapiens"))
#' orgdb <- orgdb[[1]]
#' 
#' # (4) Annotate PFAM domains
#' bp.gr <- annotatePFAM(bp.gr, orgdb)
#'
#' @export
annotatePFAM <- function(bp.gr, orgdb)
{
  up2pfam <- suppressMessages(AnnotationDbi::mapIds(orgdb, 
                                   keys = graph::nodes(bp.gr), 
                                   keytype = "UNIPROT", 
                                   column = "PFAM", 
                                   multiVals = "list"))
  graph::nodeDataDefaults(bp.gr, "PFAM") <- NA
  graph::nodeData(bp.gr, graph::nodes(bp.gr), "PFAM") <- up2pfam
  return(bp.gr)
}

# we will want to annotate these somewhat
#     ‘edgeData’: An ‘attrData’ instance for edge attributes.
#
#    ‘nodeData’: An ‘attrData’ instance for node attributes.
#
#    ‘graphData’: A ‘list’ for graph-level attributes. Only mandatory
#        list item is ‘edgemode’ which indicates whether edges are
#       ‘"directed"’ or ‘"undirected"’
#
.annotateBioplexGraph <- function(gr, df)
{
    # create maps to annotate node data:
    # 1. uniprot -> symbol
    # 2. uniprot -> entrez
    ucols <- paste0("Uniprot", c("A", "B"))
    scols <- paste0("Symbol", c("A", "B"))
    ecols <- paste0("Gene", c("A", "B"))
    
    iids <- unlist(df[,ucols])
    sids <- unlist(df[,scols])
    eids <- unlist(df[,ecols])
    
    ind <- !duplicated(iids)
    iids <- iids[ind]
    uids <- sub("-[0-9]+$", "", iids)
    up2sym <- sids[ind]
    up2eg <- eids[ind]
    names(iids) <- names(up2eg) <- names(up2sym) <- uids
    
    # graph data annotation
    # when starting off with a ordinary dfs we'll be losing the ability
    # to annotate graph-level annotation such as cell.line, version, PMID, ...
    # we might need to work with DataFrames where we have mcols and metadata
    
    # node data annotation             
    graph::nodeDataDefaults(gr, "ENTREZID") <- NA
    graph::nodeData(gr, graph::nodes(gr), "ENTREZID") <- up2eg[graph::nodes(gr)]
    graph::nodeDataDefaults(gr, "SYMBOL") <- NA
    graph::nodeData(gr, graph::nodes(gr), "SYMBOL") <- up2sym[graph::nodes(gr)]
    graph::nodeDataDefaults(gr, "ISOFORM") <- NA
    graph::nodeData(gr, graph::nodes(gr), "ISOFORM") <- iids[graph::nodes(gr)]
    
    # edge data annotation
    for(col in ucols) df[,col] <- sub("-[0-9]+$", "", df[,col])
    for(col in c("pW", "pNI", "pInt"))
    {
        graph::edgeDataDefaults(gr, col) <- numeric(0L)
        graph::edgeData(gr, df[,ucols[1]], df[,ucols[2]], col) <- df[,col]
    }
      
    return(gr)
}

#' @title Map experimental data onto a graph
#' @description Functionality for mapping experimental data stored in
#' a \code{\linkS4class{SummarizedExperiment}} onto a 
#' \code{\linkS4class{graph}} object. 
#' @param gr an object of class \code{\linkS4class{graph}}.
#' @param se an object of class \code{\linkS4class{SummarizedExperiment}}.
#' @param col.names character. Column names of \code{se} for which assay
#' data should be mapped onto the nodes of \code{gr}. Defaults to \code{NULL}
#' which will then use all column names of \code{se}.
#' @param rowdata.cols character. Column names of \code{rowData(se)} which
#' should be mapped onto the nodes of \code{gr}. Defaults to \code{NULL}
#' which will then use all column names of \code{rowData(se)}.
#' @param prefix character. Informative prefix that should be pasted together
#' with the selected \code{col.names} and \code{rowdata.cols} to allow easy
#' identification of columns of interest when mapping from multiple experimental
#' datasets.
#' @return An object of class \code{\linkS4class{graph}}. 
#' @examples
#' # (1) Obtain the latest version of the 293T PPI network ...
#' bp.293t <- getBioPlex(cell.line = "293T", version = "3.0")
#' 
#' # (2) ... and turn into a graph
#' bp.gr <- bioplex2graph(bp.293t)
#' 
#' # (3) Obtain the BioPlex3 proteome data ...
#' se <- getBioplexProteome()
#' 
#' # (4) ... and map onto the graph
#' bp.gr <- mapSummarizedExperimentOntoGraph(bp.gr, se)
#' @export
mapSummarizedExperimentOntoGraph <- function(gr, se,
                                             col.names = NULL,
                                             rowdata.cols = NULL,
                                             prefix = "")
{
    isect <- intersect(rownames(se), graph::nodes(gr))
    if(!length(isect)) stop("SummarizedExperiment (se)", " and ",
                            "graph (gr) have no node IDs in common")
    
    if(is.null(col.names)) col.names <- colnames(se)
    else if(!all(col.names %in% colnames(se)))
        stop("Invalid col.names provided")
      
    if(is.null(rowdata.cols)) rowdata.cols <- colnames(rowData(se))
    else if(!all(rowdata.cols %in% colnames(rowData(se))))
        stop("Invalid rowdata.cols provided")
    
    # map assay data onto nodes
    for(n in col.names)
    {  
        gn <- paste0(prefix, n)
        graph::nodeDataDefaults(gr, gn) <- NA
        graph::nodeData(gr, isect, gn) <- assay(se)[isect, n]
    }
    
    # map rowData onto nodes
    for(n in rowdata.cols)
    {  
        gn <- paste0(prefix, n)
        graph::nodeDataDefaults(gr, gn) <- NA
        graph::nodeData(gr, isect, gn) <- rowData(se)[isect, n]
    }
    return(gr)
}

.remapUniprotIdsBP <- function(df)
{
    stopifnot(is.data.frame(df))
    uids1 <- sub("-[0-9]+$", "", df$UniprotA)
    uids2 <- sub("-[0-9]+$", "", df$UniprotB)
    
    suppressMessages({
      ah <- AnnotationHub::AnnotationHub()
      orgdb <- AnnotationHub::query(ah, c("orgDb", "Homo sapiens"))
      orgdb <- orgdb[[1]]
    })
    
    suppressMessages({
      df$GeneA <- AnnotationDbi::mapIds(orgdb,
                                   keys = uids1,
                                   keytype = "UNIPROT",
                                   column = "ENTREZID")
      df$GeneB <- AnnotationDbi::mapIds(orgdb,
                                        keys = uids2,
                                        keytype = "UNIPROT",
                                        column = "ENTREZID")
      df$SymbolA <- AnnotationDbi::mapIds(orgdb, 
                                          keys = uids1, 
                                          keytype = "UNIPROT", 
                                          column = "SYMBOL")
      df$SymbolB <- AnnotationDbi::mapIds(orgdb, 
                                          keys = uids2, 
                                          keytype = "UNIPROT", 
                                          column = "SYMBOL")
    })
    
    return(df)
}
