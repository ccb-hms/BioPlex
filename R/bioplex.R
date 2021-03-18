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
                       remap.gene.ids = FALSE,
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

#' @export
bioplex2graph <- function(bioplex.df) 
{
    stopifnot(is.data.frame(bioplex.df))
    
    node.cols <- paste0("Uniprot", c("A", "B"))
    ftm <- as.matrix(bioplex.df[,node.cols])
    ind <- !duplicated(ftm)
    ftm <- ftm[ind,]
    bioplex.df <- bioplex.df[ind,]
    gr <- graph::ftM2graphNEL(ftm, edgemode = "directed")
    gr <- .annotateBioplexGraph(gr, bioplex.df)
    return(gr)
    
    # first approach before seeing that graph::ftM2graphNEL exists
    # maybe still useful for getting edge and node attribues onto the graph
    #
    # get protein universe (= all proteins present in the network) 
    # ie unique proteins being bait and/or prey 
    # unodes <- unique(as.vector(ftm))
    #
    # node.grid <- seq_along(unodes)
    # names(node.grid) <- unodes
    # 
    # bait <- bioplex.df[,node.cols[1]]
    # prey <- bioplex.df[,node.cols[2]]
    # 
    # indA <- node.grid[bait]
    # indB <- node.grid[prey]
    # 
    # edgeL <- split(indB, indA)
    # names(edgeL) <- unique(bait)
    # 
    # prey.only <- setdiff(prey, bait)
    # prey.only.edgeL <- replicate(length(prey.only), integer(0L), simplify = FALSE)
    # names(prey.only.edgeL) <- prey.only
    # edgeL <- c(edgeL, prey.only.edgeL)
    # graph::graphNEL(nodes = unodes, edgeL = edgeL, edgemode = "directed")
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
    
    uids <- unlist(df[,ucols])
    sids <- unlist(df[,scols])
    eids <- unlist(df[,ecols])
    
    ind <- !duplicated(uids)
    uids <- uids[ind]
    up2sym <- sids[ind]
    names(up2sym) <- uids
    up2eg <- eids[ind]
    names(up2eg) <- uids
    
    # graph data annotation
    # when starting off with a ordinary dfs we'll be losing the ability
    # to annotate graph-level annotation such as cell.line, version, PMID, ...
    # we might need to work with DataFrames where we have mcols and metadata
    
    # node data annotation             
    graph::nodeDataDefaults(gr, "ENTREZID") <- NA
    graph::nodeData(gr, graph::nodes(gr), "ENTREZID") <- up2eg[graph::nodes(gr)]
    graph::nodeDataDefaults(gr, "SYMBOL") <- NA
    graph::nodeData(gr, graph::nodes(gr), "SYMBOL") <- up2sym[graph::nodes(gr)]
    
    # edge data annotation
    for(col in c("pW", "pNI", "pInt"))
    {
        graph::edgeDataDefaults(gr, col) <- numeric(0L)
        graph::edgeData(gr, df[,ucols[1]], df[,ucols[2]], col) <- df[,col]
    }
      
    return(gr)
}
