#' merge two orgDB with keys
#' @importFrom dplyr distinct
#' @importFrom AnnotationDbi keytypes keys 
#' @importFrom RSQLite dbGetQuery dbListTables
#' @importFrom AnnotationForge makeOrgPackage
#' @importFrom utils remove.packages
#' @importFrom utils install.packages
#' @importFrom stats na.omit
#' @importFrom utils remove.packages
#' @importFrom stringr str_trim
#' @param dbleft a charater indicate the left orgDB
#' @param dbright a character indicate the right orgDB
#' @param keyleft the keytype use for merging in left orgDB
#' @param keyright the keytype use for merging in the right orgDB 
#' @param keytype the keytypes to be included in the merged orgDB ("GID","GENENAME")
#' @param keep the name of keytype you used if keyleft and keyright were not same
#' @param species the species name
#' @param author author for the annotation package
#' @param maintainer maintainer for the annotation package
#' @param tax_id taxonomy id for the species
#' @param genus genus name for the annotation package
#' @param version version number for the annotation package
#' @param pkgname package name you want to choose
#' @param install install the package or not(default: TRUE)
#' @param rebuild rebuild the package or not(default: FALSE)
#' @param outputDir temporary output path
#' @examples
#' fromKEGG(species = "hsa", anntype="KEGG")
#' fromAnnHub(species="human")
#' mergeDB("org.hsa.eg.db","org.human.eg.db",species="merge")
#' @export
#' @author Kai Guo
mergeDB<-function(dbleft,dbright,keyleft="GID",keyright="GID",keytype=NULL,keep = NULL,
                  species=NULL,author = NULL, 
                  maintainer = NULL, tax_id = NULL, genus = NULL,
                  version = NULL, pkgname=NULL,
                  install = TRUE, outputDir = NULL, rebuild = FALSE){
  ### extract keytype left
  if(!is.null(pkgname)){
    dbname <- paste0('org.', pkgname, '.eg.db')
  }else{
    dbname <- paste0('org.', species, '.eg.db')
  }
  if(isTRUE(rebuild)){
    suppressMessages(remove.packages(dbname))
  }
  if(is.null(keytype)){
    keytype=c("GID","GENENAME","SYMBOL")
  }
  #############################################
  #############################################
  keyleftl <- tolower(keyleft)
  keyrightl <- tolower(keyright)
  dbleft_name <- eval(parse(text=paste0(sub('\\.db','_dbconn()',dbleft))))
  dbright_name <- eval(parse(text=paste0(sub('\\.db','_dbconn()',dbright))))
  dblall<-dbListTables(dbleft_name)
  dbrall<-dbListTables(dbright_name)
  ##
  if(!keyleftl%in%dblall){
    keyleftl <-"genes"
  }
  if(!keyrightl%in%dbrall){
    keyrightl <-"genes"
  }
  ##
  dbl <-dbGetQuery(dbleft_name,paste0("SELECT * from"," ",keyleftl))
  dbr <-dbGetQuery(dbright_name,paste0("SELECT * from"," ",keyrightl))
  ktleft <- keytypes(eval(parse(text=dbleft)))
  ktright <- keytypes(eval(parse(text=dbright)))
  ### match the keytype in the left and right orgDB
  ksleft <- intersect(ktleft,keytype)
  ksright <- intersect(ktright,keytype)
  keys_left <- keys(eval(parse(text=dbleft)),keytype = keyleft )
  keys_right <- keys(eval(parse(text=dbleft)),keytype = keyright)
  gene2namel <- data.frame("GID" = keys_left, "GENENAME" = "")
  gene2namer <- data.frame("GID" = keys_right, "GENENAME" = "")
  gene2name <- rbind(gene2namel,gene2namer)
  if("GENENAME" %in% ksleft){
    gene2namel <- NULL
    gene2namel <- dbGetQuery(dbleft_name,"SELECT * from gene_info")
    gene2namel <- merge(dbl,gene2namel)
    gene2namel <- gene2namel[,2:3]
    colnames(gene2namel)<-c('GID','GENENAME')
  }
  if("GENENAME" %in% ksright){
    gene2namer <- NULL
    gene2namer <- dbGetQuery(dbright_name,"SELECT * from gene_info")
    gene2namer <- merge(dbr,gene2namer)
    gene2namer <- gene2namer[,2:3]
    colnames(gene2namer)<-c('GID','GENENAME')
  }
  if(!is.null(keep)){
    gene2namel[,keyleft] <- gene2namel[,1]
    gene2namer[,keyright] <- gene2namer[,1]
    colnames(gene2namel)[3] <- keep
    colnames(gene2namer)[3] <- keep
    ksleft <- setdiff(ksleft,keep)
    ksright <- setdiff(ksright,keep)
  }else{
    gene2namel[,keyleft] <- gene2namel[,1]
    gene2namer[,keyright] <- gene2namer[,1]
  }
  geneinfo<-rbind(gene2namel,gene2namer)
  geneinfo<-na.omit(geneinfo)
  geneinfo<-distinct(geneinfo)
  colnames(geneinfo)[1:2]<-c('GID','GENENAME')
  geneinfo$GID<-str_trim(geneinfo$GID,side = "both")
  geneinfo$GENENAME<-str_trim(geneinfo$GENENAME,side = "both")
  if(nrow(geneinfo)>1){
    geneinfo <- geneinfo[geneinfo$GENENAME!="",]
  }
  gene2gol <- data.frame("GID" = geneinfo$GID[1], "GO" = "", 
                        "EVIDENCE" = "IEA")
  gene2gor <- data.frame("GID" = geneinfo$GID[1], "GO" = "", 
                         "EVIDENCE" = "IEA")
  if("GO" %in% ksleft){
    gene2gol <- NULL
    ####eval(parse(text=paste0("dbListTables(org.mac.eg","_dbconn())")))
    gene2gol <- dbGetQuery(dbleft_name,"SELECT * from go_all")
    gene2gol <- merge(dbl,gene2gol)
    gene2gol <- gene2gol[,2:4]
  }
  if("GO" %in% ksright){
    gene2gor <- NULL
    gene2gor <- dbGetQuery(dbright_name,"SELECT * from go_all")
    gene2gor <- merge(dbr,gene2gor)
    gene2gor <- gene2gor[,2:4]
  }
  gene2go<-rbind(gene2gol,gene2gor)
  gene2go<-na.omit(gene2go)
  gene2go<-distinct(gene2go)
  colnames(gene2go)<-c('GID','GO','EVIDENCE')
  if(nrow(gene2go)>1){
    gene2go<-gene2go[gene2go$GO!="",]
  }
  gene2pathl <- data.frame("GID" = geneinfo$GID[1],"PATH" = "01100")
  gene2pathr <- data.frame("GID" = geneinfo$GID[1],"PATH" = "01100")
  if("PATH" %in% ksleft){
    gene2pathl <- NULL
    if("path" %in% dblall){
      gene2pathl <- dbGetQuery(dbleft_name,"SELECT * from path")
    }else{
      gene2pathl <- dbGetQuery(dbleft_name,"SELECT * from kegg")
    }
    gene2pathl <- merge(dbl,gene2pathl)
    gene2pathl <- gene2pathl[,2:3]
  }
  if("PATH" %in% ksright){
    gene2pathr <- NULL
    if("path" %in% dbrall){
      gene2pathr <- dbGetQuery(dbright_name,"SELECT * from path")
    }else{
      gene2pathr <- dbGetQuery(dbright_name,"SELECT * from kegg")
    }
    gene2pathr <- merge(dbl,gene2pathr)
    gene2pathr <- gene2pathr[,2:3]
    
  }
  gene2path<-rbind(gene2pathl,gene2pathr)
  gene2path<-na.omit(gene2path)
  gene2path<-distinct(gene2path)
  colnames(gene2path)<-c('GID','PATH')
  if(nrow(gene2path)>1){
    gene2path<-gene2path[gene2path$PATH!="",]
  }
  gene2kol <- data.frame("GID" = geneinfo$GID[1],"KO" = "")
  gene2kor <- data.frame("GID" = geneinfo$GID[1],"KO" = "")
  if("KO" %in% ksleft){
    gene2kol <- NULL
    gene2kol <- dbGetQuery(dbleft_name,"SELECT * from ko")
    gene2kol <- merge(dbl,gene2kol)
    gene2kol <- gene2kol[,2:3]
  }
  if("KO" %in% ksright){
    gene2kor <- NULL
    gene2kor <- dbGetQuery(dbright_name,"SELECT * from ko")
    gene2kor <- merge(dbl,gene2kor)
    gene2kor <- gene2kor[,2:3]
  }
  gene2ko<-rbind(gene2kol,gene2kor)
  gene2ko<-na.omit(gene2ko)
  gene2ko<-distinct(gene2ko)
  colnames(gene2ko)<-c('GID','KO')
  if(nrow(gene2path)>1){
    gene2ko<-gene2ko[gene2ko$KO!="",]
  }
  gene2refseql <- data.frame("GID" = geneinfo$GID[1], "REFSEQ" = "")
  gene2refseqr <- data.frame("GID" = geneinfo$GID[1], "REFSEQ" = "")
  if("REFSEQ" %in% ksleft){
    gene2refseql <- NULL
    gene2refseql <- dbGetQuery(dbleft_name,"SELECT * from refseq")
    gene2refseql <- merge(dbl,gene2refseql)
    gene2refseql <- gene2refseql[,2:3]
  }
  if("REFSEQ" %in% ksright){
    gene2refseqr <- NULL
    gene2refseqr <- dbGetQuery(dbright_name,"SELECT * from refseq")
    gene2refseqr <- merge(dbl,gene2refseqr)
    gene2refseqr <- gene2refseqr[,2:3]
  }
  gene2refseq<-rbind(gene2refseql,gene2refseqr)
  colnames(gene2refseq)<-c('GID','REFSEQ')
  if(nrow(gene2refseq)>1){
    gene2refseq<-gene2refseq[gene2refseq$REFSEQ!="",]
  }
  if(keep == "REFSEQ"){
    gene2refseq <- geneinfo[,c(1,3)]
    geneinfo <- geneinfo[,1:2]
  }
  gene2refseq<-na.omit(gene2refseq)
  gene2refseq<-distinct(gene2refseq)
  gene2symboll <- data.frame("GID" = geneinfo$GID[1], "SYMBOL" = "")
  gene2symbolr <- data.frame("GID" = geneinfo$GID[1], "SYMBOL" = "")
  if("SYMBOL" %in% ksleft){
    gene2symboll <- NULL
    gene2symboll <- dbGetQuery(dbleft_name,"SELECT * from symbol")
    gene2symboll <- merge(dbl,gene2symboll)
    gene2symboll <- gene2symboll[,2:3]
  }
  if("SYMBOL" %in% ksright){
    gene2symbolr <- NULL
    gene2symbolr <- dbGetQuery(dbright_name,"SELECT * from symbol")
    gene2symbolr <- merge(dbl,gene2symbolr)
    gene2symbolr <- gene2symbolr[,2:3]
  }
  gene2symbol<-rbind(gene2symboll,gene2symbolr)
  colnames(gene2symbol)<-c('GID','SYMBOL')
  if(nrow(gene2symbol)>1){
    gene2symbol<-gene2symbol[gene2symbol$SYMBOL!="",]
  }
  if(keep == "SYMBOL"){
    gene2symbol <- geneinfo[,c(1,3)]
    geneinfo <- geneinfo[,1:2]
  }
  gene2symbol<-na.omit(gene2symbol)
  gene2symbol<-distinct(gene2symbol)
  #
  gene2ensembll <- data.frame("GID" = geneinfo$GID[1], "ENSEMBL" = "")
  gene2ensemblr <- data.frame("GID" = geneinfo$GID[1], "ENSEMBL" = "")
  if("ENSEMBL" %in% ksleft){
    gene2ensembll <- NULL
    gene2ensembll <- dbGetQuery(dbleft_name,"SELECT * from ensembl")
    gene2ensembll <- merge(dbl,gene2ensembll)
    gene2ensembll <- gene2ensembll[,2:3]
  }
  if("ENSEMBL" %in% ksright){
    gene2ensemblr <- NULL
    gene2ensemblr <- dbGetQuery(dbright_name,"SELECT * from ensembl")
    gene2ensemblr <- merge(dbl,gene2ensemblr)
    gene2ensemblr <- gene2ensemblr[,2:3]
  }
  gene2ensembl<-rbind(gene2ensembll,gene2ensemblr)
  colnames(gene2ensembl)<-c('GID','ENSEMBL')
  if(nrow(gene2ensembl)>1){
    gene2ensembl<-gene2ensembl[gene2ensembl$ENSEMBL!="",]
  }
  #
  if(keep == "ENSEMBL"){
    gene2ensembl <- geneinfo[,c(1,3)]
    geneinfo <- geneinfo[,1:2]
  }
  gene2ensembl<-na.omit(gene2ensembl)
  gene2ensembl<-distinct(gene2ensembl)
  #
  gene2entrezidl <- data.frame("GID" = geneinfo$GID[1], "ENTREZID" = "")
  gene2entrezidr <- data.frame("GID" = geneinfo$GID[1], "ENTREZID" = "")
  if("ENTREZID" %in% ksleft){
    gene2entrezidl <- NULL
    gene2entrezidl <- dbGetQuery(dbleft_name,"SELECT * from entrezid")
    gene2entrezidl <- merge(dbl,gene2entrezidl)
    gene2entrezidl <- gene2entrezidl[,2:3]
  }
  if("ENTREZID" %in% ksright){
    gene2entrezidr <- NULL
    gene2entrezidr <- dbGetQuery(dbright_name,"SELECT * from entrezid")
    gene2entrezidr <- merge(dbl,gene2entrezidr)
    gene2entrezidr <- gene2entrezidr[,2:3]
  }
  gene2entrezid<-rbind(gene2entrezidl,gene2entrezidr)
  colnames(gene2entrezid)<-c('GID','ENTREZID')
  if(nrow(gene2entrezid)>1){
    gene2entrezid<-gene2entrezid[gene2entrezid$ENTREZID!="",]
  }
  #
  if(keep == "ENTREZID"){
    gene2entrezid <- geneinfo[,c(1,3)]
    geneinfo <- geneinfo[,1:2]
  }
  gene2entrezid<-na.omit(gene2entrezid)
  gene2entrezid<-distinct(gene2entrezid)
  gene2pfaml <- data.frame("GID" = geneinfo$GID[1], "PFAM" = "")
  gene2pfamr <- data.frame("GID" = geneinfo$GID[1], "PFAM" = "")
  if("PFAM" %in% ksleft){
    gene2pfaml <- NULL
    gene2pfaml <- dbGetQuery(dbleft_name,"SELECT * from pfam")
    gene2pfaml <- merge(dbl,gene2pfaml)
    gene2pfaml <- gene2pfaml[,2:3]
  }
  if("PFAM" %in% ksright){
    gene2pfamr <- NULL
    gene2pfamr <- dbGetQuery(dbright_name,"SELECT * from pfam")
    gene2pfamr <- merge(dbl,gene2pfamr)
    gene2pfamr <- gene2pfamr[,2:3]
  }
  gene2pfam<-rbind(gene2pfaml,gene2pfamr)
  gene2pfam<-na.omit(gene2pfam)
  gene2pfam<-distinct(gene2pfam)
  colnames(gene2pfam)<-c('GID','PFAM')
  if(nrow(gene2pfam)>1){
    gene2pfam<-gene2pfam[gene2pfam$PFAM!="",]
  }
  gene2interprol <- data.frame("GID" = geneinfo$GID[1], "INTERPRO" = "")
  gene2interpror <- data.frame("GID" = geneinfo$GID[1], "INTERPRO" = "")
  if("INTERPRO" %in% ksleft){
    gene2interprol <- NULL
    gene2interprol <- dbGetQuery(dbleft_name,"SELECT * from interpro")
    gene2interprol <- merge(dbl,gene2interprol)
    gene2interprol <- gene2interprol[,2:3]
  }
  if("INTERPRO" %in% ksright){
    gene2interpror <- NULL
    gene2interpror <- dbGetQuery(dbright_name,"SELECT * from interpro")
    gene2interpror <- merge(dbl,gene2interpror)
    gene2interpror <- gene2interpror[,2:3]
  }
  gene2interpro<-rbind(gene2interprol,gene2interpror)
  gene2interpro<-na.omit(gene2interpro)
  gene2interpro<-distinct(gene2interpro)
  colnames(gene2interpro)<-c('GID','INTERPRO')
  if(nrow(gene2interpro)>1){
    gene2interpro<-gene2interpro[gene2interpro$INTERPRO!="",]
  }
  #
  gene2reactl <- data.frame("GID" = geneinfo$GID[1], "REACTOME" = "")
  gene2reactr <- data.frame("GID" = geneinfo$GID[1], "REACTOME" = "")
  if("REACTOME" %in% ksleft){
    gene2reactl <- NULL
    gene2reactl <- dbGetQuery(dbleft_name,"SELECT * from reactome")
    gene2reactl <- merge(dbl,gene2reactl)
    gene2reactl <- gene2reactl[,2:3]
  }
  if("REACTOME" %in% ksright){
    gene2reactr <- NULL
    gene2reactr <- dbGetQuery(dbright_name,"SELECT * from reactome")
    gene2reactr <- merge(dbl,gene2reactr)
    gene2reactr <- gene2reactr[,2:3]
  }
  gene2react<-rbind(gene2reactl,gene2reactr)
  gene2react<-na.omit(gene2react)
  gene2react<-distinct(gene2react)
  colnames(gene2react)<-c('GID','REACTOME')
  if(nrow(gene2react)>1){
    gene2react<-gene2react[gene2react$REACTOME!="",]
  }
  #
  gene2biocycl <- data.frame("GID" = geneinfo$GID[1], "BIOCYC" = "")
  gene2biocycr <- data.frame("GID" = geneinfo$GID[1], "BIOCYC" = "")
  if("BIOCYC" %in% ksleft){
    gene2biocycl <- NULL
    gene2biocycl <- dbGetQuery(dbleft_name,"SELECT * from biocyc")
    gene2biocycl <- merge(dbl,gene2biocycl)
    gene2biocycl <- gene2biocycl[,2:3]
  }
  if("BIOCYC" %in% ksright){
    gene2biocycr <- NULL
    gene2biocycr <- dbGetQuery(dbright_name,"SELECT * from biocyc")
    gene2biocycr <- merge(dbl,gene2biocycr)
    gene2biocycr <- gene2biocycr[,2:3]
  }
  gene2biocyc<-rbind(gene2biocycl,gene2biocycr)
  gene2biocyc<-na.omit(gene2biocyc)
  gene2biocyc<-distinct(gene2biocyc)
  colnames(gene2biocyc)<-c('GID','BIOCYC')
  if(nrow(gene2biocyc)>1){
    gene2biocyc<-gene2biocyc[gene2biocyc$BIOCYC!="",]
  }
  ###
  gene2kdl <- data.frame("GID" = geneinfo$GID[1], "KEGGDISEASE" = "")
  gene2kdr <- data.frame("GID" = geneinfo$GID[1], "KEGGDISEASE" = "")
  if("KEGGDISEASE" %in% ksleft){
    gene2kdl <- NULL
    gene2kdl <- dbGetQuery(dbleft_name,"SELECT * from disease")
    gene2kdl <- merge(dbl,gene2kdl)
    gene2kdl <- gene2kdl[,2:3]
  }
  if("KEGGDISEASE" %in% ksright){
    gene2kdr <- NULL
    gene2kdr <- dbGetQuery(dbright_name,"SELECT * from disease")
    gene2kdr <- merge(dbl,gene2kdr)
    gene2kdr <- gene2kdr[,2:3]
  }
  gene2kd<-rbind(gene2kdl,gene2kdr)
  gene2kd<-na.omit(gene2kd)
  gene2kd<-distinct(gene2kd)
  colnames(gene2kd)<-c('GID','KEGGDISEASE')
  if(nrow(gene2kd)>1){
    gene2kd<-gene2kd[gene2kd$KEGGDISEASE!="",]
  }
  ####
  gene2gadl <- data.frame("GID" = geneinfo$GID[1], "GAD" = "")
  gene2gadr <- data.frame("GID" = geneinfo$GID[1], "GAD" = "")
  if("GAD" %in% ksleft){
    gene2gadl <- NULL
    gene2gadl <- dbGetQuery(dbleft_name,"SELECT * from gad")
    gene2gadl <- merge(dbl,gene2gadl)
    gene2gadl <- gene2gadl[,2:3]
  }
  if("GAD" %in% ksright){
    gene2gadr <- NULL
    gene2gadr <- dbGetQuery(dbright_name,"SELECT * from gad")
    gene2gadr <- merge(dbl,gene2gadr)
    gene2gadr <- gene2gadr[,2:3]
  }
  gene2gad<-rbind(gene2gadl,gene2gadr)
  gene2gad<-na.omit(gene2gad)
  gene2gad<-distinct(gene2gad)
  colnames(gene2gad)<-c('GID','GAD')
  if(nrow(gene2gad)>1){
    gene2gad<-gene2gad[gene2gad$GAD!="",]
  }
  ###
  gene2fundol <- data.frame("GID" = geneinfo$GID[1], "FUNDO" = "")
  gene2fundor <- data.frame("GID" = geneinfo$GID[1], "FUNDO" = "")
  if("FUNDO" %in% ksleft){
    gene2fundol <- NULL
    gene2fundol <- dbGetQuery(dbleft_name,"SELECT * from fundo")
    gene2fundol <- merge(dbl,gene2fundol)
    gene2fundol <- gene2fundol[,2:3]
  }
  if("FUNDO" %in% ksright){
    gene2fundor <- NULL
    gene2fundor <- dbGetQuery(dbright_name,"SELECT * from fundo")
    gene2fundor <- merge(dbl,gene2fundor)
    gene2fundor <- gene2fundor[,2:3]
  }
  gene2fundo<-rbind(gene2fundol,gene2fundor)
  gene2fundo<-na.omit(gene2fundo)
  gene2fundo<-distinct(gene2fundo)
  colnames(gene2fundo)<-c('GID','FUNDO')
  if(nrow(gene2fundo)>1){
    gene2fundo<-gene2fundo[gene2fundo$FUNDO!="",]
  }
  ###
  if(is.null(author)){
    author <-"myself"
  }
  if(is.null(maintainer)){
    maintainer <- "mysel<myself@gmail.com>"
  }
  if(is.null(tax_id)){
    tax_id <- "123"
  }
  if(is.null(species)){
    species <- species
  }
  if(!is.null(pkgname)){
    species <- pkgname
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
  geneinfo<-na.omit(geneinfo)
  geneinfo<-distinct(geneinfo)
  package <- suppressWarnings(makeOrgPackage(gene_info = geneinfo,
                                             symbol = gene2symbol,
                                             entrezid = gene2entrezid,
                                             refseq = gene2refseq,
                                             ensembl = gene2ensembl,
                                             go = gene2go,
                                             path = gene2path,
                                             ko = gene2ko,
                                             pfam = gene2pfam,
                                             interpro = gene2interpro,
                                             reactome = gene2react,
                                             biocyc = gene2biocyc,
                                             disease = gene2kd,
                                             gad = gene2gad,
                                             fundo = gene2fundo,
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
