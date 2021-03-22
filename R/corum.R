############################################################
# 
# author: Ludwig Geistlinger
# date: 2021-02-16 09:20:02
# 
# descr: download CORUM data 
# 
############################################################

#' @title Obtain CORUM protein complex data
#' @description Functionality for retrieving the CORUM protein complex data.
#' Available complex collections include:
#' \itemize{
#' \item complete set of complexes,
#' \item core set of complexes,
#' \item complexes with splice variants.
#' } See references.
#' @param set character. Valid options include:
#' \itemize{ 
#' \item \code{"all"}: complete set of complexes,
#' \item \code{"core"}: core set of complexes,
#' \item \code{"splice"}: complexes with splice variants.
#' }
#' Defaults to \code{"all"}.
#' @param organism character. Use \code{NULL} to not subset by organism.
#' Defaults to \code{"Human"} which restricts the data to human protein complexes
#' only.
#' @param remap.uniprot.ids logical. Should the protein-to-gene mappings from CORUM
#' (i.e. UNIPROT-to-SYMBOL and UNIPROT-to-ENTREZID) be updated using Bioc annotation
#' functionality? Currently only supported in combination with \code{organism = "Human"}.
#' Defaults to \code{FALSE} which will then keep the mappings provided by CORUM.
#' See details.
#' @param cache logical. Should a locally cached version used if available?
#' Defaults to \code{TRUE}.
#' @return A \code{data.frame}. 
#' @references CORUM: \url{http://mips.helmholtz-muenchen.de/corum/#download}
#' @examples
#' # Obtain the core set of CORUM complexes
#' core <- getCorum(set = "core")
#' @importFrom utils download.file relist unzip
#' @export
getCorum <- function(set = c("all", "core", "splice"),
                     organism = "Human",
                     remap.uniprot.ids = FALSE,
                     cache = TRUE)
{
    corum.url <- "http://mips.helmholtz-muenchen.de/corum/download"
    set <- match.arg(set)
    rname <- paste("corum", 
                   set,
                   ifelse(is.null(organism), "all", organism),
                   ifelse(remap.uniprot.ids, "remapped", "original"),
                   sep = "-")

    # should a cache version be used?
    if(cache)
    {
        corum <- .getResourceFromCache(rname)
        if(!is.null(corum)) return(corum)        
    }

    # download, unzip, and read in
    set <- paste0(set, "Complexes.txt.zip")
    download.file(file.path(corum.url, set), destfile = set)
    unzip(set)
    set <- sub(".zip$", "", set)
    corum <- read.delim(set)
    
    # organism
    if(!is.null(organism))
    {
        stopifnot(organism %in% corum$Organism)
        corum <- subset(corum, Organism %in% organism)
    }
    
    # remap gene ids
    if(remap.uniprot.ids) corum <- .remapUniprotIds(corum)
    
    # clean up & cache
    .cacheResource(corum, rname)
    file.remove(c(set, paste0(set, ".zip")))
    return(corum) 
}

.remapUniprotIds <- function(df)
{
    if(any(df$Organism != "Human")) 
        stop("Gene ID re-mapping is currently only supported for human")
    if(!requireNamespace("AnnotationHub"))
        stop("Please install the 'AnnotationHub' package for gene ID mapping")

    # relevant column names
    id.col <- "subunits.UniProt.IDs."
    eg.col <- "subunits.Entrez.IDs."
    sym.col <- "subunits.Gene.name."
    
    # obtain latest ensembl mappings from AnnotationHub
    suppressMessages({
        ah <- AnnotationHub::AnnotationHub()
        ahdb <- AnnotationHub::query(ah, c("orgDb", "Homo sapiens"))
        orgdb <- ahdb[[length(ahdb)]]
    })
    
    # map Uniprot IDs to Entrez IDs and to symbols
    unip <- strsplit(df[[id.col]], " ?; ?")
    ulnip <- unlist(unip)
    uunip <- unique(ulnip)
    
    suppressMessages({
        egs <- AnnotationDbi::mapIds(orgdb,
                                      keys = uunip,
                                      keytype = "UNIPROT",
                                      column = "ENTREZID")
        syms <- AnnotationDbi::mapIds(orgdb, 
                                      keys = uunip, 
                                      keytype = "UNIPROT", 
                                      column = "SYMBOL")
    })
    
    # replace corresponding columns in corum df
    egl <- relist(unname(egs[ulnip]), unip)
    df[[eg.col]] <- vapply(egl, paste, character(1), collapse = ";") 
    syl <- relist(unname(syms[ulnip]), unip)
    df[[sym.col]] <- vapply(syl, paste, character(1), collapse = ";")
    
    return(df)
}


#' @export
corum2list <- function(corum.df, 
                       subunit.id.type = c("UNIPROT", "ENTREZID"))
{
    stopifnot(is.data.frame(corum.df))
    subunit.id.type <- match.arg(subunit.id.type)
    subunit.id.type <- ifelse(subunit.id.type == "UNIPROT", "UniProt", "Entrez")
    id.col <- paste0("subunits.", subunit.id.type, ".IDs.")
    
    gs <- strsplit(corum.df[[id.col]], " ?; ?")
    complex.name <- gsub(" ", "_", corum.df$ComplexName)
    names(gs) <- paste0("CORUM", corum.df$ComplexID, "_", complex.name)
    return(gs)
}

#' basic data structure for sets of complexes or sets of pull-downs:
#' a list of graph instances
#' @export
corum2graphlist <- function(corum.df, 
                            subunit.id.type = c("UNIPROT", "ENTREZID"))
{
    cl <- corum2list(corum.df, subunit.id.type)
    gl <- lapply(cl, .completelyConnected)
    gl <- .annotateGraphList(gl, corum.df)
    return(gl)
}

# @param glist a list of graphs
# @param a gene 
# @param SYMBOL and returns a logical vector indicating
# @return which graphs have a node with the input gene SYMBOL
#' @export
hasSubunit <- function(glist, subunit, id.type = "SYMBOL") 
{
    .hasX <- function(x) grep(subunit, 
                               unlist(graph::nodeData(x, graph::nodes(x), id.type)), 
                               ignore.case = TRUE)
    res <- lapply(glist, .hasX) 
    lengths(res) > 0
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
.annotateGraphList <- function(glist, df)
{
    # FIXME: Someone needs to fix the graph class accessors for graph attributes
    #- we should not need to access the slot directly
    go.ids <- strsplit(df[["GO.ID"]], " ?; ?")
    entrez.ids <- strsplit(df[["subunits.Entrez.IDs."]], " ?; ?")
    symbols <- strsplit(df[["subunits.Gene.name."]], " ?; ?")
    for(i in seq_along(glist))
    {
        current <- glist[[i]]
        gd <- list(ComplexID = df[i, "ComplexID"],
                   ComplexName = df[i, "ComplexName"], 
                   GO.ID = go.ids[[i]], 
                   PubMed.ID = df[i, "PubMed.ID"])
        glist[[i]]@graphData <- c(current@graphData, gd)
        graph::nodeDataDefaults(glist[[i]], "ENTREZID") <- NA
        
        if(length(entrez.ids[[i]]) == graph::numNodes(current))
            graph::nodeData(glist[[i]], 
                            graph::nodes(glist[[i]]), 
                            "ENTREZID") <- entrez.ids[[i]]
        graph::nodeDataDefaults(glist[[i]], "SYMBOL") <- NA
        if(length(symbols[[i]]) == graph::numNodes(current))
            graph::nodeData(glist[[i]], 
                            graph::nodes(glist[[i]]), 
                            "SYMBOL") <- symbols[[i]]
        
    }
    return(glist)
}

# since our complexes are completely connected we will need a little
# helper function to create the instances
# @param nodes character vector of unique node names 
# @return instance of the graphNEL class where all nodes are connected 
# to all other nodes with undirected edges
.completelyConnected <- function(nodes) 
{
    if(!is.character(nodes)) stop("require character vector for nodes")
    if(length(nodes) == 0 ) stop("require at least one node")
    if(length(nodes) == 1 ) 
        return(graph::graphNEL(nodes = nodes, edgeL = list()))
    
    edL <- vector("list", length = length(nodes))
    names(edL) <- nodes
    eInds <- seq_along(nodes)
    for(i in eInds)
        edL[[i]] <- list(edges = eInds[-i])
    return(graph::graphNEL(nodes = nodes, edgeL = edL))
}
