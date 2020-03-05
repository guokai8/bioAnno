#' build Own annotation database with user defined annotation file
#' @importFrom AnnotationForge makeOrgPackage
#' @importFrom utils install.packages
#' @importFrom stats na.omit
#' @param geneinfo gene information table with two columns
#'                as default("GID","DESCRIPTION")
#' @param gene2go Gene Onotoly information for  genes
#' @param gene2path KEGG Pathway or KO information for genes
#' @param author author for the annotation package
#' @param maintainer maintainer for the annotation package
#' @param tax_id taxonomy id for the species
#' @param genus genus for the species
#' @param version version for the annotation package
#' @param species species name(common name,kegg.species.code or scientifc name)
#' @param install install the package or not(default: TRUE)
#' @param outputDir temporary output path
#' @export
#' @examples
#' require(bioAnno)
#' data(ath)
#' fromOwn(geneinfo = ath)
#' @return annotation package
#' @author Kai Guo
fromOwn <- function(geneinfo = geneinfo, gene2go = NULL, gene2path = NULL,
        version = NULL, maintainer = NULL, author = NULL, outputDir = NULL,
        tax_id = NULL, genus = NULL, species = NULL, install = TRUE){

    cat("Please make sure you have Gene Ontology and KEGG pathway
        or KO data.frame ready.\n")

    if(is.null(geneinfo)){
        stop("You must have Gene information table")
    }
    if(colnames(geneinfo)[1] != "GID"){
    stop("The first column name must be GID")
    }
    geneinfo <- geneinfo[!duplicated(geneinfo), ]
    geneinfo <- na.omit(geneinfo)
    if(!is.null(gene2go)){
        if(ncol(gene2go) == 2){
        gene2go$EVIDENCE <- "IEA"
        colnames(gene2go)[c(1,2)] <- c("GID", "GO")
    }else{
        colnames(gene2go) <- c("GID", "GO", "EVIDENCE")
    }
        gene2go <- gene2go[!duplicated(gene2go), ]
        gene2go <- na.omit(gene2go)
    }else{
    gene2go <- data.frame("GID" = geneinfo$GID,
        "GO" = "GO:0008150", "EVIDENCE" = "IEA")
    }
    if(!is.null(gene2path)){
        if(ncol(gene2path) != 2){
        stop("Pathway datafram must have only two columns")
    }
    colnames(gene2path) <- c("GID", "PATH")
    gene2path <- gene2path[!duplicated(gene2path), ]
    gene2path <- na.omit(gene2path)
    }
    if(is.null(version)){
    version <- "0.0.1"
    }
    if(is.null(tax_id)){
    tax_id <- "xxx"
    }
    if(is.null(author)){
    author <- "myself"
    }
    if(is.null(maintainer)){
    maintainer <- "myself<myself@email.com>"
    }
    if(is.null(genus)){
    genus <- ""
    }
    if(is.null(species)){
    species <- "species"
    }
    if(is.null(outputDir)){
    outputDir <- tempdir()
    }
    package <- makeOrgPackage(
    gene_info = geneinfo,
    go = gene2go,
    path = gene2path,
    version = version,
    maintainer = maintainer,
    author = author,
    outputDir = outputDir,
    tax_id = tax_id,
    genus = genus,
    species = species,
    goTable = "go"
    )
    if(isTRUE(install)){
    install.packages(package, repos = NULL, type = "source")
    unlink(package, recursive = TRUE)
    }
}
