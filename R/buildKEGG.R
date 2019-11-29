#' make annotation database from KEGG and GO from NCBI
#' @importFrom KEGGREST keggLink
#' @importFrom KEGGREST keggList
#' @importFrom AnnotationForge makeOrgPackage
#' @param species species name(common name,kegg.species.code or scientifc name)
#' @param author author for the annotation package
#' @param maintainer maintainer
#' @param tax_id taxonomy id for the species
#' @param genus genus for the species
#' @param version version
#' @param install install the package or not(default: TRUE)
#' @author Kai Guo
#' @export
#'
fromKEGG<-function(species="ath",anntype=NULL,buildall=TRUE,author=NULL,
                   maintainer=NULL,tax_id=NULL,genus=NULL,version=NULL,
                   install=TRUE,outputDir=NULL,rebuild=FALSE){
  cat("##############################################################################################\n")
  cat("The bioAnno package downloads and uses KEGG data. Non-academic uses may require a KEGG license\nagreement (details at http://www.kegg.jp/kegg/legal.html)\n")
  cat("The Gene Ontology are downloaded from NCBI.\n")
  cat("##############################################################################################\n")
  dbinfo <- .get.species.info(species)
  species <- dbinfo["kegg.code"]
  dbname <- paste0('org.',species,'.eg.db')
  if (require(dbname,character.only=TRUE)&& !isTRUE(rebuild)){
    suppressMessages(require(dbname,character.only = T,quietly = T))
  }else{
  geneinfo <- data.frame()
  gene2path <- data.frame()
  gene2ko <- data.frame()
  tmp<-keggList(species)
  geneinfo <- data.frame("GID"=sub(paste0(species,":"),'',names(tmp)),"GENAME"= tmp)
  gene2go <- .extratGO(taxid=dbinfo["tax.id"])
  gene2go <- gene2go[!duplicated(gene2go),]
  tmp <- keggLink('pathway', species)
  gene2path <- data.frame("GID"=sub(paste0(species,":"),'',names(tmp)),
                          "PATH"= sub(species,'',sub('path:','',tmp)))
  tmp <- keggLink('ko', species)
  gene2ko <- data.frame("GID"= sub(paste0(species,":"),'',names(tmp)),"KO"= sub('ko:','',tmp))

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
    outputDir<-tempdir()
  }
  package<-makeOrgPackage(
    gene_info=geneinfo,
    path=gene2path,
    ko=gene2ko,
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
    install.packages(package,repos = NULL,type="source")
  }
  }
}


