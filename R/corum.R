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
                   #ifelse(is.null(organism), "all", organism),
                   #ifelse(remap.uniprot.ids, "remapped", "original"),
                   sep = "-")

    # should a cache version be used?
    if(cache) set.file <- .getResourceFromCache2(rname)
    if(!cache || is.null(set.file))        
    {
        # download, unzip, and read in
        set.file <- paste0(set, "Complexes.txt.zip")
        set.url <- file.path(corum.url, set.file)
        set.file <- .cacheResource2(rname, set.url)#, download = FALSE)
    }

    set.file <- unzip(set.file, exdir = dirname(set.file))
    corum <- read.delim(set.file)
    file.remove(set.file)
    
    # organism
    if(!is.null(organism))
    {
        stopifnot(organism %in% corum$Organism)
        corum <- subset(corum, Organism %in% organism)
    }
    
    # remap gene ids
    if(remap.uniprot.ids) corum <- .remapUniprotIds(corum)
    
    return(corum) 
}

#' @title Represent CORUM protein complex data as a simple list
#' @description Functionality for storing CORUM protein complex data in a
#' \code{list}.
#' @param corum.df A \code{data.frame} storing the CORUM protein complex data.
#' Typically obtained via \code{\link{getCorum}}.
#' @param subunit.id.type character. Supported options include \code{"UNIPROT"}
#' (default) and \code{"ENTREZID"}.
#' @return A \code{list} with an entry for each complex. Each entry is a 
#' character vector of subunit IDs. 
#' @references CORUM: \url{http://mips.helmholtz-muenchen.de/corum/#download}
#' @examples
#'  # (1) Obtain the core set of CORUM complexes ...
#'  core <- getCorum(set = "core")
#'  
#'  # (2) ... turn into a list
#'  core.list <- corum2list(core)
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

#' @title Represent CORUM protein complex data as a list of graph instances
#' @description Functionality for storing CORUM protein complex data in a
#' \code{list} of \code{graph} instances.
#' @param corum.df A \code{data.frame} storing the CORUM protein complex data.
#' Typically obtained via \code{\link{getCorum}}.
#' @param subunit.id.type character. Supported options include \code{"UNIPROT"}
#' (default) and \code{"ENTREZID"}.
#' @return A \code{list} with an entry for each complex. Each entry is an 
#' object of class \code{graphNEL} connecting all subunit IDs with each other
#' by undirected edges. 
#' @references CORUM: \url{http://mips.helmholtz-muenchen.de/corum/#download}
#' @examples
#'  # (1) Obtain the core set of CORUM complexes ...
#'  core <- getCorum(set = "core")
#'  
#'  # (2) ... turn into a list of graphs
#'  core.glist <- corum2graphlist(core)
#' @export
corum2graphlist <- function(corum.df, 
                            subunit.id.type = c("UNIPROT", "ENTREZID"))
{
    cl <- corum2list(corum.df, subunit.id.type)
    gl <- lapply(cl, .completelyConnected)
    gl <- .annotateGraphList(gl, corum.df)
    return(gl)
}

#' @title Identify CORUM complexes that have a subunit of interest
#' @description Screens a \code{list} of \code{graph} instances storing 
#' CORUM protein complex data for a subunit of choice.
#' @param glist A \code{list} of \code{graph}s storing CORUM complexes.
#' Typically obtained via \code{\link{corum2graphlist}}. 
#' @param subunit character. A gene ID corresponding to the subunit of interest. 
#' @param id.type character. Gene ID type of the given subunit. Defaults to \code{"SYMBOL"}. 
#' @return A logical vector indicating which graphs have a node with the given 
#' subunit.
#' @examples
#' # (1) Obtain the core set of CORUM complexes ...
#' core <- getCorum(set = "core")
#' 
#' # (2) ... turn into a list of graphs ...
#' core.glist <- corum2graphlist(core)
#' 
#' # (3) .. check for a particular subunit of interest
#' has.cdk2 <- hasSubunit(core.glist, subunit = "CDK2")
#' @export
hasSubunit <- function(glist, subunit, id.type = "SYMBOL") 
{
    .hasX <- function(x) grep(subunit, 
                               unlist(graph::nodeData(x, graph::nodes(x), id.type)), 
                               ignore.case = TRUE)
    res <- lapply(glist, .hasX) 
    lengths(res) > 0
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
        # AnnotationHub::query(ah, c("orgDb", "Homo sapiens"))
        orgdb <- ah[["AH92581"]]
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
