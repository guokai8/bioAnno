#'make annotation database from NCBI
#' @importFrom AnnotationForge makeOrgPackage
#' @param species species name
#' @param author author for the annotation package
#' @param maintainer maintainer
#' @param tax_id taxonomy id for the species
#' @param genus genus for the species
#' @param version version
#' @param plant plant or animal species (TRUE/FALSE)
#' @param install install the package or not(default: TRUE)
#' @author Kai Guo
#' @export
fromNCBI <- function(species="ath",author=NULL,
                               maintainer=NULL,tax_id=NULL,genus=NULL,version=NULL,
                               install=TRUE,outputDir=NULL){
  dbinfo <- .get.species.info(species=species)
  dbname <- paste0('org.',species,'.eg.db')
  if (require(dbname,character.only=TRUE)){
    suppressMessages(require(dbname,character.only = T,quietly = T))
  }else{
  geneinfo <- .extratGene(taxid = dbinfo['tax.id'])
  gene2symbol<-geneinfo[,c("GID","SYMBOL")]
  gene2symbol[!duplicated(gene2symbol),]
  geneinfo <- geneinfo[,c("GID","DESCRIPTION")]
  geneinfo <- geneinfo[!duplicated(geneinfo),]
  gene2go <- .extratGO(taxid = dbinfo['tax.id'])
  if(is.null(version)){
    version="0.0.1"
  }
  if(is.null(tax_id)){
    tax_id=dbinfo["tax.id"]
  }
  if(is.null(author)){
    author="myself"
  }
  if(is.null(maintainer)){
    maintainer="myself<myself@email.com>"
  }
  if(is.null(genus)){
    genus=dbinfo["scientific.name"]
  }
  if(is.null(species)){
    species=species
  }
  if(is.null(outputDir)){
    outputDir="."
  }
  package<-makeOrgPackage(
    gene_info=geneinfo,
    gene2symbol=gene2symbol,
    go=gene2go,
    maintainer=maintainer,
    author=author,
    outputDir = outputDir,
    tax_id=tax_id,
    genus="",
    species=species,
    version=version,
    goTable="go"
  )
  tmp <- NULL
  if(isTRUE(install)){
    package=sub('.*\\/','',package)
    install.packages(package,repos = NULL,type="source")
  }
}
}
