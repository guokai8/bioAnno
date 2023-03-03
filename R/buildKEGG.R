#' @title make annotation database from KEGG and GO from NCBI
#' @importFrom KEGGREST keggLink
#' @importFrom KEGGREST keggList
#' @importFrom AnnotationForge makeOrgPackage
#' @importFrom utils remove.packages
#' @importFrom utils install.packages
#' @importFrom stats na.omit
#' @importFrom utils data
#' @param species species name(common name,kegg.species.code or scientifc name)
#' @param anntype the type of function annotation(GO,KEGG,PFAM,InterPro)
#'             you want get from ensemble
#' @param author author for the annotation package
#' @param maintainer maintainer for the annotation package
#' @param tax_id taxonomy id for the species
#' @param genus genus for the species
#' @param version version for the annotation package
#' @param install install the package or not(default: TRUE)
#' @param rebuild rebuild the package or not(default: FALSE)
#' @param outputDir temporary output path
#' @examples
#' fromKEGG(species = "eco", install = FALSE)
#' @author Kai Guo
#' @return annotation package
#' @export
#'
fromKEGG <- function(species="ath", anntype=c("KEGG","GO"), author=NULL,
                maintainer=NULL,tax_id=NULL,genus=NULL,version=NULL,
                install=TRUE,outputDir=NULL,rebuild=FALSE){
    cat("#########################################################################\n")
    cat("The bioAnno package downloads and uses KEGG data.Non-academic uses may
require a KEGG license agreement (details at http://www.kegg.jp/kegg/legal.html)\n")
    cat("The Gene Ontology are downloaded from NCBI.\n")
    cat("#########################################################################\n")
    dbinfo <- .get.species.info(species)
    species <- dbinfo["kegg.code"]
    dbname <- paste0('org.', species, '.eg.db')
    if(isTRUE(rebuild)){
        suppressMessages(remove.packages(dbname))
    }
    if(is_installed(dbname)){
        suppressMessages(requireNamespace(dbname,quietly = TRUE))
        cat("You alreay had the annotation package: ", dbname, " \n")
    }else{
#  if (require(dbname,character.only=TRUE) & !isTRUE(rebuild)){
#    suppressMessages(require(dbname,character.only = T,quietly = T))
#  }else{
    geneinfo <- data.frame()
    gene2path <- data.frame()
    gene2ko <- data.frame()
    tmp <- keggList(species)
    geneinfo <- data.frame("GID" = sub(paste0(species,":"), '',
        names(tmp)), "GENENAME" = tmp)
    rownames(geneinfo) <- NULL
    if("GO" %in% anntype){
      gene2go <- .extratGO(taxid = dbinfo["tax.id"])
      gene2go <- gene2go[!duplicated(gene2go), ]
    }
    if(nrow(gene2go) == 0){
        cat("No Gene Ontology information available !\n")
        gene2go <- data.frame("GID" = geneinfo$GID,
        "GO" = "GO:0008150", "EVIDENCE" = "IEA")
    }
    if(species=="ath"){
        if(!exists("ath")) data(ath)
        gene2go$GID <- ath[gene2go$GID, 1]
        gene2go <- na.omit(gene2go)
    }
    tmp <- keggLink('pathway', species)
    gene2path <- data.frame("GID" = sub(paste0(species,":"), '', names(tmp)),
        "PATH"= sub(species,'',sub('path:','',tmp)))
    tmp <- keggLink('ko', species)
    gene2ko <- data.frame("GID" = sub(paste0(species, ":"), '',
        names(tmp)), "KO" = sub('ko:', '', tmp))

    if(is.null(version)){
        version <- "0.0.1"
    }
    if(is.null(tax_id)){
        tax_id <- dbinfo["tax.id"]
    }
    if(is.null(author)){
        author <- "myself"
    }
    if(is.null(maintainer)){
        maintainer <- "myself<myself@email.com>"
    }
    if(is.null(genus)){
        genus <- dbinfo["scientific.name"]
    }
    if(is.null(species)){
        species <- species
    }
    if(is.null(outputDir)){
        outputDir <- tempdir()
    }
    package <- suppressWarnings(makeOrgPackage(
    gene_info = geneinfo,
    path = gene2path,
    ko = gene2ko,
    go = gene2go,
    maintainer = maintainer,
    author = author,
    outputDir = outputDir,
    tax_id = tax_id,
    genus = "",
    species = species,
    version = version,
    verbose = FALSE,
    goTable = "go"
    ))
    tmp <- NULL
    if(isTRUE(install)){
        install.packages(package, repos = NULL, type="source")
        unlink(package, recursive = TRUE)
    }else{
        .show.path(package)
        .show.tables(package)
        return(package)
    }
    }
}


