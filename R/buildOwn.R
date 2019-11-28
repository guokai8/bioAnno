#' build Own annotation database with user defined annotation file
#' @importFrom AnnotationForge makeOrgPackage
#' @param geneinfo gene information table with two columns as default("GID","DESCRIPTION")
#' @param entrezid gene information table with corresponding ENTREZID id
#' @param refseq gene information table with corresponding REFSEQ id
#' @param gene2go Gene Onotoly information for each gene
#' @param author author for the annotation package
#' @param maintainer maintainer
#' @param tax_id taxonomy id for the species
#' @param genus genus for the species
#' @param version version
#' @param install install the package or not(default: TRUE)
fromOwn<-function(geneinfo=geneinfo,entrezid=NULL,refseq=NULL,gene2go=NULL,gene2path=NULL,pfam=NULL,interpro=NULL,reactome=NULL,
                   version=NULL,maintainer=NULL,author=NULL, outputDir = NULL,tax_id=NULL,genus=NULL,species=NULL){

  if(is.null(geneinfo)){
    stop("You must have Gene information table")
  }
  if(colnames(geneinfo)[1]!="GID"){
    stop("The first column name must be GID")
  }
  geneinfo<-geneinfo[!duplicated(geneinfo),]
  geneinfo<-na.omit(geneinfo)
  if(!is.null(entrezid)){
    if(ncol(entrezid)!=2){
      stop("entrezid dataframe must have only had two columns")
    }
    colnames(entrezid)<-c("GID","ENTREZID")
    entrezid<-entrezid[!duplicated(entrezid),]
    entrezid<-na.omit(entrezid)
    entrezid$ENTREZID<-as.character(entrezid$ENTREZID)
  }
  if(!is.null(refseq)){
    if(ncol(refseq)!=2){
      stop("refseq dataframe must have only had two columns")
    }
    colnames(refseq)<-c("GID","REFSEQ")
    refseq<-refseq[!duplicated(refseq),]
    refseq<-na.omit(refseq)
  }
  if(!is.null(gene2go)){
    if(ncol(gene2go)==2){
      gene2go$EVIDENCE<-"IEA"
      colnames(gene2go)[1:2]<-c("GID","GO")
    }else{
      colnames(gene2go)<-c("GID","GO","EVIDENCE")
    }
    gene2go<-gene2go[!duplicated(gene2go),]
    gene2go<-na.omit(gene2go)
  }else{
    gene2go <- data.frame("GID"=geneinfo$GID,"GO"="GO:0008150","EVIDENCE"="IEA")
  }
  if(!is.null(gene2path)){
    if(ncol(gene2path)!=2){
      stop("Pathway datafram must have only two columns")
    }
    colnames(gene2path)<-c("GID","PATH")
    gene2path<-gene2path[!duplicated(gene2path),]
    gene2path<-na.omit(gene2path)
  }
  if(!is.null(pfam)){
    if(ncol(pfam)!=2){
      stop("PFAM datafram must have only two columns")
    }
    colnames(pfam)<-c("GID","PFAM")
    pfam<-pfam[!duplicated(pfam),]
    pfam<-na.omit(pfam)
  }
  if(!is.null(interpro)){
    if(ncol(interpro)!=2){
      stop("InterPro datafram must have only two columns")
    }
    colnames(interpro)<-c("GID","INTERPRO")
    interpro<-interpro[!duplicated(interpro),]
    interpro<-na.omit(interpro)
  }
  if(!is.null(reactome)){
    if(ncol(reactome)!=2){
      stop("Reactome datafram must have only two columns")
    }
    colnames(reactome)<-c("GID","REACTOME")
    reactome<-reactome[!duplicated(reactome),]
    reactome<-na.omit(reactome)
  }
  if(is.null(version)){
    version="0.0.1"
  }
  if(is.null(tax_id)){
    tax_id="xxx"
  }
  if(is.null(author)){
    author="myself"
  }
  if(is.null(maintainer)){
    maintainer="myself<myself@email.com>"
  }
  if(is.null(genus)){
    genus=""
  }
  if(is.null(species)){
    species="species"
  }
  if(is.null(outputDir)){
    outputDir="."
  }
  package<-makeOrgPackage(
  gene_info=geneinfo,
  entrezid=gene2entrezid,
  refseq=gene2ref,
  go=gene2go,
  path=gene2path,
  pfam=pfam,
  interpro=interpro,
  reactome=reactome,
  version=version,
  maintainer=maintainer,
  author=author,
  outputDir = outputDir,
  tax_id=tax_id,
  genus=genus,
  species=species,
  goTable="go"
  )
  if(isTRUE(install)){
    package <- sub('.*\\/','',package)
    install.packages(package,repos = NULL,type="source")
  }
}
