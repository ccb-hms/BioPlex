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
#' @param cache logical. Should a locally cached version used if available?
#' Defaults to \code{TRUE}.
#' @return A \code{data.frame}. 
#' @references CORUM: \url{http://mips.helmholtz-muenchen.de/corum/#download}
#' @examples
#' # Obtain the core set of CORUM complexes
#' core <- getCorum(set = "core")
#' @importFrom utils download.file unzip
#' @export
getCorum <- function(set = c("all", "core", "splice"),
                     cache = TRUE)
{
    corum.url <- "http://mips.helmholtz-muenchen.de/corum/download"
    set <- match.arg(set)
    rname <- paste("corum", set, sep = "-")

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

    # clean up & cache
    .cacheResource(corum, rname)
    file.remove(c(set, paste0(set, ".zip")))
    return(corum) 
}

#' @export
corum2list <- function(corum.df, 
                       organism = "Human",
                       subunit.id.type = c("UNIPROT", "ENTREZID"))
{
    subunit.id.type <- match.arg(subunit.id.type)
    stopifnot(organism %in% corum.df$Organism)
    
    corum.df <- subset(corum.df, Organism == organism)
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
                            organism = "Human",
                            subunit.id.type = c("UNIPROT", "ENTREZID"))
{
    cl <- corum2list(corum.df, organism, subunit.id.type)
    gl <- lapply(cl, .completelyConnected)
    gl <- .annotateGraphList(gl, corum.df)
    return(gl)
}

# @param glist a list of graphs
# @param a gene 
# @param SYMBOL and returns a logical vector indicating
# @return which graphs have a node with the input gene SYMBOL
.hasSubunit <- function(glist, subunit, id.type = "SYMBOL") 
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
.annotateGraphList <- function(glist, df, remap.node.ids = FALSE)
{
    # FIXME: Someone needs to fix the graph class accessors for graph attributes
    #- we should not need to access the slot directly
    go.ids <- strsplit(df[["GO.ID"]], " ?; ?")
    for(i in seq_along(glist)) 
        glist[[i]]@graphData <- c(glist[[i]]@graphData, 
                                           ComplexID = df[i, "ComplexID"],
                                           ComplexName = df[i, "ComplexName"], 
                                           GO.ID = go.ids[[i]], 
                                           PubMed.ID = df[i, "PubMed.ID"])
    if(remap.node.ids) .annotateNodeIds(glist)
    return(glist)
}

.annotateNodeIds <- function(glist, annotation, from = "UNIPROT", to = "SYMBOL")
{
    for(i in seq_along(glist))
    {
        graph::nodeDataDefaults(glist[[i]], to) <- NA 
        #graph::nodeData(glist[[i]], graph::nodes(glist[[i]]), to) <- toEG[[i]]
    }    
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