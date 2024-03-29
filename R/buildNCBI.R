#' @title build annotation database from NCBI
#' @importFrom AnnotationForge makeOrgPackage
#' @importFrom utils remove.packages
#' @param species species name
#' @param author author for the annotation package
#' @param maintainer maintainer for the annotation package
#' @param tax_id taxonomy id for the species
#' @param genus genus for the species
#' @param version version for the annotation package
#' @param install install the package or not(default: TRUE)
#' @param pkgname package name you want to choose
#' @param rebuild rebuild the package or not(default: FALSE)
#' @param outputDir temporary output path
#' @examples
#' \donttest{
#' ## build annoataion package for Ecoli
#' fromNCBI(species = "eco", install = FALSE)
#' }
#' @author Kai Guo
#' @return annotation package
#' @export
fromNCBI <- function(species = "ath", author = NULL,
        maintainer = NULL, tax_id = NULL, genus=NULL, version = NULL,
        install = TRUE, pkgname = NULL, outputDir=NULL, rebuild = FALSE){
    dbinfo <- .get.species.info(species=species)
    species <- dbinfo["kegg.code"]
    dbname <- paste0('org.',species,'.eg.db')
    if(isTRUE(rebuild)){
        suppressMessages(remove.packages(dbname))
    }
    if(is_installed(dbname)){
        suppressMessages(requireNamespace(dbname,quietly = TRUE))
        cat("You alreay had the annotation package: ",dbname," \n")
    }else{
    #  if (require(dbname,character.only=TRUE)){
    #    suppressMessages(require(dbname,character.only = T,quietly = T))
    #  }else{
        geneinfo <- .extratGene(taxid = dbinfo['tax.id'])
        gene2entrezid <- data.frame("GID"=geneinfo$GID,"ENTREZID"= geneinfo$GID)
        gene2symbol<-geneinfo[,c("GID","SYMBOL")]
        gene2symbol[!duplicated(gene2symbol),]
        geneinfo <- geneinfo[,c("GID","DESCRIPTION")]
        geneinfo <- geneinfo[!duplicated(geneinfo),]
        gene2go <- .extratGO(taxid = dbinfo['tax.id'])
    if(nrow(gene2go)==0){
        cat("No Gene Ontology information available !\n")
        gene2go <- data.frame("GID" = geneinfo$GID, "GO" = "GO:0008150", 
                    "EVIDENCE" = "IEA")
    }
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
    if(!is.null(pkgname)){
        species <- pkgname
    }
    if(is.null(outputDir)){
        outputDir <- tempdir()
    }
    package <- suppressWarnings(makeOrgPackage(
    gene_info = geneinfo,
    gene2symbol = gene2symbol,
    entrezid = gene2entrezid,
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
        install.packages(package, repos = NULL, type = "source")
        unlink(package,recursive = TRUE)
    }else{
        .show.path(package)
        .show.tables(package)
        #return(package)
    }
    }
}
