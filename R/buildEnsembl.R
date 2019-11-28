#' build annotation from ensembl
#' @title build annotation from ensembl
#' @importFrom biomaRt useMart
#' @importFrom biomaRt useDataset
#' @importFrom biomaRt getBM
#' @importFrom AnnotationForge makeOrgPackage
#' @param host the ensemble API host,for plant you can use plants.ensembl.org and for human and other species you can use uswest.ensembl.org
#' @param species the sepcies you want to search, you can use showplant to get the species name
#' @param anntype the type of function annotation(GO,KEGG,PFAM,InterPro) you want get from ensemble
#' @param buildall include all prossbile annoation type listed in Ensembl
#' @param author author for the annotation package
#' @param maintainer maintainer
#' @param tax_id taxonomy id for the species
#' @param genus genus name
#' @param version version number(xx.xx.xx)
#' @param plant plant or animal species (TRUE/FALSE)
#' @param install install the package or not(default: TRUE)
#' @author Kai Guo
#' @export
fromEnsembl<-function(species="Arabidopsis t",host="uswest.ensembl.org",
                       anntype=NULL,buildall=TRUE,author=NULL,
                       maintainer=NULL,tax_id=NULL,genus=NULL,version=NULL,plant=FALSE,
                       install=TRUE,outputDir=NULL){
  if(isTRUE(plant)){
     host="plants.ensembl.org"
     mart=useMart("plants_mart",host=host)
  }else{
     mart=useMart("ENSEMBL_MART_ENSEMBL",host=host)
  }
  dbinfo<-.getmartdb(species,mart)
 # dbname1 <- paste0('org.',strsplit(species," ")[[1]][1],'.eg.db')
 #  if (require(dbname1,character.only=TRUE)){
 #    suppressMessages(require(dbname1,character.only = T,quietly = T))
#  }else{
  dbname=as.character(dbinfo$dbname)
  dataset<-useDataset(dbname,mart=mart)
  if(is.null(anntype)){
    if(isTRUE(buildall)){
      anntype=c("GO","KEGG","Reactome","PFAM","InterPro")
    }else{
      stop("You need to specific anntation type!\n")
    }
  }
  geneinfo<-getBM(attributes = c("ensembl_gene_id","description"),filters ="chromosome_name",values = as.vector(dbinfo$chr_info$name),dataset)
 # geneinfo<-geneinfo[nchar(geneinfo[,2])>1,]
  colnames(geneinfo)<-c("GID","GENENAME")
  geneinfo<-na.omit(geneinfo)
  gene2entrezid<-getBM(attributes = c("ensembl_gene_id","entrezgene_id"),filters ="chromosome_name",values = as.vector(dbinfo$chr_info$name),dataset)
  gene2entrezid<-gene2entrezid[nchar(gene2entrezid[,2])>1,]
  colnames(gene2entrezid)<-c("GID","ENTREZID")
  gene2entrezid<-na.omit(gene2entrezid)
  gene2entrezid[,2]<-as.character(gene2entrezid[,2])
  gene2ref<-getBM(attributes = c("ensembl_gene_id","refseq_dna"),filters ="chromosome_name",values = as.vector(dbinfo$chr_info$name),dataset)
  gene2ref<-gene2ref[nchar(gene2ref[,2])>1,]
  colnames(gene2ref)<-c("GID","REFSEQ")
  gene2ref<-na.omit(gene2ref)
  if("GO"%in%anntype){
    gene2go<-getBM(attributes = c("ensembl_gene_id","go_id","go_linkage_type"),filters ="chromosome_name",values = as.vector(dbinfo$chr_info$name),dataset)
    gene2go<-gene2go[nchar(gene2go[,2])>1,]
    gene2go<-gene2go[nchar(gene2go[,3])>1,]
    colnames(gene2go)<-c("GID","GO","EVIDENCE")
    gene2go<-na.omit(gene2go)
    gene2go<-gene2go[!duplicated(gene2go),]
  }
  if("KEGG"%in%anntype){
    gene2path<-getBM(attributes = c("ensembl_gene_id","kegg_enzyme"),filters ="chromosome_name",values = as.vector(dbinfo$chr_info$name),dataset)
    gene2path[,2]<-sub('\\+.*','',gene2path[,2])
    gene2path<-gene2path[nchar(gene2path[,2])>1,]
    colnames(gene2path)<-c("GID","PATH")
    gene2path<-na.omit(gene2path)
    gene2path<-gene2path[!duplicated(gene2path),]
  }
  if("PFAM"%in%anntype){
    gene2pfam<-getBM(attributes = c("ensembl_gene_id","pfam"),filters ="chromosome_name",values = as.vector(dbinfo$chr_info$name),dataset)
    gene2pfam<-gene2pfam[nchar(gene2pfam[,2])>1,]
    colnames(gene2pfam)<-c("GID","PFAM")
    gene2pfam<-na.omit(gene2pfam)
  }
  if("InterPro"%in%anntype){
    gene2interpro<-getBM(attributes = c("ensembl_gene_id","interpro"),filters ="chromosome_name",values = as.vector(dbinfo$chr_info$name),dataset)
    gene2interpro<-gene2interpro[nchar(gene2interpro[,2])>1,]
    colnames(gene2interpro)<-c("GID","INTERPRO")
    gene2interpro<-na.omit(gene2interpro)
  }
  if("Reactome"%in%anntype){
    if(isTRUE(plant)){
      gene2reactome<-getBM(attributes = c("ensembl_gene_id","plant_reactome_pathway"),filters ="chromosome_name",values = as.vector(dbinfo$chr_info$name),dataset)
    }else{
      gene2reactome<-getBM(attributes = c("ensembl_gene_id","reactome"),filters ="chromosome_name",values = as.vector(dbinfo$chr_info$name),dataset)
    }
    gene2reactome<-gene2reactome[nchar(gene2reactome[,2])>1,]
    colnames(gene2reactome)<-c("GID","REACTOME")
    gene2reactome<-na.omit(gene2reactome)
  }
  if(is.null(author)){
    author="myself"
  }
  if(is.null(maintainer)){
    maintainer="mysel<myself@gmail.com>"
  }
  if(is.null(tax_id)){
    tax_id="123"
  }
  if(is.null(version)){
    version="0.0.1"
  }
  if(is.null(genus)){
    genus=""
  }
  if(is.null(outputDir)){
    outputDir<-"."
  }
  species=gsub(' .*','',species)
  package<-makeOrgPackage(gene_info=geneinfo,
                 entrezid=gene2entrezid,
                 refseq=gene2ref,
                 go=gene2go,
                 path=gene2path,
                 pfam=gene2pfam,
                 interpro=gene2interpro,
                 reactome=gene2reactome,
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
    package=sub('.*\\/','',package)
    install.packages(package,repos = NULL,type="source")
 # }
  }
}
##'
##' @importFrom  biomaRt useMart
##' @importFrom  biomaRt listDatasets
##' @param host Ensembl host site
##' @author Kai Guo
##' @export
listSpecies<-function(host="uswest.ensembl.org",plant=FALSE){
  cat("You could choose different host to get high speed!\n")
  if(isTRUE(plant)){
    host="plants.ensembl.org"
    mart=useMart("plants_mart",host=host)
  }else{
    cat("Ensembl US West: uswest.ensembl.org\nEnsembl US East: useast.ensembl.org\nEnsembl Asia: asia.ensembl.org\n" )
    mart=useMart("ENSEMBL_MART_ENSEMBL",host=host)
  }
  res<-listDatasets(mart)
  colnames(res)[2]<-"species"
  res
}

