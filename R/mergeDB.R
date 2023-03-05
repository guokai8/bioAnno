#' merge two orgDB with keys
#' @importFrom dplyr distinct
#' @importFrom AnnotationDbi keytypes keys 
#' @importFrom AnnotationForge makeOrgPackage
#' @importFrom utils remove.packages
#' @importFrom utils install.packages
#' @importFrom stats na.omit
#' @importFrom utils remove.packages
#' @param dbleft the left orgDB
#' @param dbright the right orgDB
#' @param keyleft the keytype use for merging in left orgDB
#' @param keyright the keytype use for merging in the right orgDB 
#' @param keytype the keytypes to be included in the merged orgDB ("GID","GENENAME")
#' @param species the species name
#' @param author author for the annotation package
#' @param maintainer maintainer for the annotation package
#' @param tax_id taxonomy id for the species
#' @param genus genus name for the annotation package
#' @param version version number for the annotation package
#' @param install install the package or not(default: TRUE)
#' @param rebuild rebuild the package or not(default: FALSE)
#' @param outputDir temporary output path
#' @examples
#' fromKEGG(species = "hsa", anntype="KEGG")
#' fromAnnHub(species="human")
#' mergeDB(org.hsa.eg.db,org.human.eg.db,species="merge")
#' @export
#' @author Kai Guo
mergeDB<-function(dbleft,dbright,keyleft="GID",keyright="GID",keytype=NULL,
                  species=NULL,author = NULL,
                  maintainer = NULL, tax_id = NULL, genus = NULL,
                  version = NULL, 
                  install = TRUE, outputDir = NULL, rebuild = FALSE){
  ### extract keytype left
  dbname <- paste0('org.',species,'.eg.db')
  if(isTRUE(rebuild)){
    suppressMessages(remove.packages(dbname))
  }
  if(is.null(keytype)){
    keytype=c("GID","GENENAME","SYMBOL")
  }
  keytype <-setdiff(keytype,c(keyleft,keyright))
  ktleft <- keytypes(dbleft)
  ktright <- keytypes(dbright)
  ### match the keytype in the left and right orgDB
  ksleft <- intersect(ktleft,keytype)
  ksright <- intersect(ktright,keytype)
  keys_left <- keys(dbright,keytype=keyleft)
  keys_right <- keys(dbright,keytype=keyright)
  gene2namel <- data.frame("GID" = keys_left, "GENENAME" = "")
  gene2namer <- data.frame("GID" = keys_right, "GENENAME" = "")
  gene2name <- data.frame("GID" = keys_right, "GENENAME" = "")
  if("GENENAME" %in% ksleft){
    gene2namel <- NULL
    gene2namel <- AnnotationDbi::select(dbleft,keys=keys_left,keytype=keyleft,columns=c("GENENAME"))
  }
  if("GENENAME" %in% ksright){
    gene2namer <- NULL
    gene2namer <- AnnotationDbi::select(dbright,keys=keys_right,keytype=keyright,columns=c("GENENAME"))
  }
  geneinfo<-rbind(gene2namel,gene2namer)
  geneinfo<-na.omit(geneinfo)
  geneinfo<-distinct(geneinfo)
  colnames(geneinfo)<-c('GID','GENENAME')
  gene2gol <- data.frame("GID" = geneinfo$GID[1], "GO" = "GO:0008150", 
                        "EVIDENCE" = "IEA")
  gene2gor <- data.frame("GID" = geneinfo$GID[1], "GO" = "GO:0008150", 
                         "EVIDENCE" = "IEA")
  if("GO" %in% ksleft){
    gene2gol <- NULL
    gene2gol <- AnnotationDbi::select(dbleft,keys=keys_left,keytype=keyleft,columns=c("GOALL","EVIDENCE"))
  }
  if("GO" %in% ksright){
    gene2gor <- NULL
    gene2gor <- AnnotationDbi::select(dbright,keys=keys_right,keytype=keyright,columns=c("GOALL","EVIDENCE"))
  }
  gene2go<-rbind(gene2gol,gene2gor)
  gene2go<-na.omit(gene2go)
  gene2go<-distinct(gene2go)
  colnames(gene2go)<-c('GID','GO','EVIDENCE')
  gene2pathl <- data.frame("GID" = geneinfo$GID[1],"PATH" = "01100")
  gene2pathr <- data.frame("GID" = geneinfo$GID[1],"PATH" = "01100")
  if("PATH" %in% ksleft){
    gene2pathl <- NULL
    gene2pathl <- AnnotationDbi::select(dbleft,keys=keys_left,keytype=keyleft,columns=c("PATH"))
  }
  if("PATH" %in% ksright){
    gene2pathr <- NULL
    gene2pathr <- AnnotationDbi::select(dbright,keys=keys_right,keytype=keyright,columns=c("PATH"))
  }
  gene2path<-rbind(gene2pathl,gene2pathr)
  gene2path<-na.omit(gene2path)
  gene2path<-distinct(gene2path)
  colnames(gene2path)<-c('GID','PATH')
  gene2kol <- data.frame("GID" = geneinfo$GID[1],"KO" = "")
  gene2kor <- data.frame("GID" = geneinfo$GID[1],"KO" = "")
  if("KO" %in% ksleft){
    gene2kol <- NULL
    gene2kol <- AnnotationDbi::select(dbleft,keys=keys_left,keytype=keyleft,columns=c("KO"))
  }
  if("KO" %in% ksright){
    gene2kor <- NULL
    gene2kor <- AnnotationDbi::select(dbright,keys=keys_right,keytype=keyright,columns=c("KO"))
  }
  gene2ko<-rbind(gene2kol,gene2kor)
  gene2ko<-na.omit(gene2ko)
  gene2ko<-distinct(gene2ko)
  colnames(gene2ko)<-c('GID','KO')
  gene2refseql <- data.frame("GID" = geneinfo$GID[1], "REFSEQ" = "")
  gene2refseqr <- data.frame("GID" = geneinfo$GID[1], "REFSEQ" = "")
  if("REFSEQ" %in% ksleft){
    gene2refseql <- NULL
    gene2refseql <- AnnotationDbi::select(dbleft,keys=keys_left,keytype=keyleft,columns=c("REFSEQ"))
  }
  if("REFSEQ" %in% ksright){
    gene2refseqr <- NULL
    gene2refseqr <- AnnotationDbi::select(dbright,keys=keys_right,keytype=keyright,columns=c("REFSEQ"))
  }
  gene2refseq<-rbind(gene2refseql,gene2refseqr)
  gene2refseq<-na.omit(gene2refseq)
  gene2refseq<-distinct(gene2refseq)
  colnames(gene2refseq)<-c('GID','REFSEQ')
  gene2symboll <- data.frame("GID" = geneinfo$GID[1], "SYMBOL" = "")
  gene2symbolr <- data.frame("GID" = geneinfo$GID[1], "SYMBOL" = "")
  if("SYMBOL" %in% ksleft){
    gene2symboll <- NULL
    gene2symboll <- AnnotationDbi::select(dbleft,keys=keys_left,keytype=keyleft,columns=c("SYMBOL"))
  }
  if("SYMBOL" %in% ksright){
    gene2symbolr <- NULL
    gene2symbolr <- AnnotationDbi::select(dbright,keys=keys_right,keytype=keyright,columns=c("SYMBOL"))
  }
  gene2symbol<-rbind(gene2symboll,gene2symbolr)
  gene2symbol<-na.omit(gene2symbol)
  gene2symbol<-distinct(gene2symbol)
  colnames(gene2symbol)<-c('GID','SYMBOL')
  #
  gene2ensembll <- data.frame("GID" = geneinfo$GID[1], "ENSEMBL" = "")
  gene2ensemblr <- data.frame("GID" = geneinfo$GID[1], "ENSEMBL" = "")
  if("ENSEMBL" %in% ksleft){
    gene2ensembll <- NULL
    gene2ensembll <- AnnotationDbi::select(dbleft,keys=keys_left,keytype=keyleft,columns=c("ENSEMBL"))
  }
  if("ENSEMBL" %in% ksright){
    gene2ensemblr <- NULL
    gene2ensemblr <- AnnotationDbi::select(dbright,keys=keys_right,keytype=keyright,columns=c("ENSEMBL"))
  }
  gene2ensembl<-rbind(gene2ensembll,gene2ensemblr)
  gene2ensembl<-na.omit(gene2ensembl)
  gene2ensembl<-distinct(gene2ensembl)
  colnames(gene2ensembl)<-c('GID','ENSEMBL')
  #
  gene2entrezidl <- data.frame("GID" = geneinfo$GID[1], "ENTREZID" = "")
  gene2entrezidr <- data.frame("GID" = geneinfo$GID[1], "ENTREZID" = "")
  if("ENTREZID" %in% ksleft){
    gene2entrezidl <- NULL
    gene2entrezidl <- AnnotationDbi::select(dbleft,keys=keys_left,keytype=keyleft,columns=c("ENTREZID"))
  }
  if("ENTREZID" %in% ksright){
    gene2entrezidr <- NULL
    gene2entrezidr <- AnnotationDbi::select(dbright,keys=keys_right,keytype=keyright,columns=c("ENTREZID"))
  }
  gene2entrezid<-rbind(gene2entrezidl,gene2entrezidr)
  gene2entrezid<-na.omit(gene2entrezid)
  gene2entrezid<-distinct(gene2entrezid)
  colnames(gene2entrezid)<-c('GID','ENTREZID')
  gene2pfaml <- data.frame("GID" = geneinfo$GID[1], "PFAM" = "")
  gene2pfamr <- data.frame("GID" = geneinfo$GID[1], "PFAM" = "")
  if("PFAM" %in% ksleft){
    gene2pfaml <- NULL
    gene2pfaml <- AnnotationDbi::select(dbleft,keys=keys_left,keytype=keyleft,columns=c("PFAM"))
  }
  if("PFAM" %in% ksright){
    gene2pfamr <- NULL
    gene2pfamr <- AnnotationDbi::select(dbright,keys=keys_right,keytype=keyright,columns=c("PFAM"))
  }
  gene2pfam<-rbind(gene2pfaml,gene2pfamr)
  gene2pfam<-na.omit(gene2pfam)
  gene2pfam<-distinct(gene2pfam)
  colnames(gene2pfam)<-c('GID','PFAM')
  gene2interprol <- data.frame("GID" = geneinfo$GID[1], "INTERPRO" = "")
  gene2interpror <- data.frame("GID" = geneinfo$GID[1], "INTERPRO" = "")
  if("INTERPRO" %in% ksleft){
    gene2interprol <- NULL
    gene2interprol <- AnnotationDbi::select(dbleft,keys=keys_left,keytype=keyleft,columns=c("INTERPRO"))
  }
  if("INTERPRO" %in% ksright){
    gene2interpror <- NULL
    gene2interpror <- AnnotationDbi::select(dbright,keys=keys_right,keytype=keyright,columns=c("INTERPRO"))
  }
  gene2interpro<-rbind(gene2interprol,gene2interpror)
  gene2interpro<-na.omit(gene2interpro)
  gene2interpro<-distinct(gene2interpro)
  colnames(gene2interpro)<-c('GID','INTERPRO')
  gene2reactl <- data.frame("GID" = geneinfo$GID[1], "REACTOME" = "")
  gene2reactr <- data.frame("GID" = geneinfo$GID[1], "REACTOME" = "")
  if("REACTOME" %in% ksleft){
    gene2reactl <- NULL
    gene2reactl <- AnnotationDbi::select(dbleft,keys=keys_left,keytype=keyleft,columns=c("REACTOME"))
  }
  if("REACTOME" %in% ksright){
    gene2reactr <- NULL
    gene2reactr <- AnnotationDbi::select(dbright,keys=keys_right,keytype=keyright,columns=c("REACTOME"))
  }
  gene2react<-rbind(gene2reactl,gene2reactr)
  gene2react<-na.omit(gene2react)
  gene2react<-distinct(gene2react)
  colnames(gene2react)<-c('GID','REACTOME')
  if(is.null(author)){
    author <-"myself"
  }
  if(is.null(maintainer)){
    maintainer <- "mysel<myself@gmail.com>"
  }
  if(is.null(tax_id)){
    tax_id <- "123"
  }
  if(is.null(version)){
    version <- "0.0.1"
  }
  if(is.null(genus)){
    genus <- ""
  }
  if(is.null(outputDir)){
    outputDir <- tempdir()
  }
  species <- gsub(' .*', '', species)
  package <- suppressWarnings(makeOrgPackage(gene_info = geneinfo,
                                             symbol = gene2symbol,
                                             entrezid = gene2entrezid,
                                             refseq = gene2refseq,
                                             go = gene2go,
                                             path = gene2path,
                                             ko = gene2ko,
                                             pfam = gene2pfam,
                                             interpro = gene2interpro,
                                             reactome = gene2react,
                                             ensembl = gene2ensembl,
                                             version = version,
                                             maintainer = maintainer,
                                             author = author,
                                             outputDir = outputDir,
                                             tax_id = tax_id,
                                             genus = genus,
                                             species = species,
                                             verbose = FALSE,
                                             goTable = "go"
  ))
  if(isTRUE(install)){
    install.packages(package, repos = NULL, type = "source")
    unlink(package, recursive = TRUE)
  }else{
    .show.path(package)
    .show.tables(package)
    #return(package)
  }
}
