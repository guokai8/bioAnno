#' build Own annotation database with user defined annotation file
#' @importFrom AnnotationForge makeOrgPackage
#' @importFrom utils install.packages
#' @importFrom stats na.omit
#' @param geneinfo gene information table with two columns
#'                as default("GID","DESCRIPTION")
#' @param gene2go Gene Onotoly information for  genes
#' @param gene2path KEGG Pathway information for genes
#' @param gene2symbol SYMBOL information for genes
#' @param gene2refseq REFSEQ or KO information for genes
#' @param gene2ensembl ENSEMBL or KO information for genes
#' @param gene2pfam PFAM information for genes
#' @param gene2reactome REACTOME Pathway or KO information for genes
#' @param gene2ko KO information for genes
#' @param gene2interpro INTERPRO information for genes
#' @param gene2entrezid ENTREZID information for genes
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
#' ## build your own annotation for Arabidopsis thaliana
#' data(ath)
#' fromOwn(geneinfo = ath, install = FALSE)
#' @return annotation package
#' @author Kai Guo
fromOwn <- function(geneinfo = geneinfo, gene2go = NULL, gene2path = NULL, 
                    gene2symbol = NULL, gene2refseq = NULL,  gene2ensembl = NULL,
                    gene2pfam = NULL, gene2reactome= NULL, gene2ko = NULL,
                    gene2interpro = NULL, gene2entrezid= NULL,
        version = NULL, maintainer = NULL, author = NULL, outputDir = NULL,
        tax_id = NULL, genus = NULL, species = NULL, install = TRUE){

    cat("Please make sure you have Gene Ontology and KEGG pathway
        or KO data.frame ready.\n")

    if(is.null(geneinfo)){
        stop("You must have Gene information table")
    }
    colnames(geneinfo)[1] <- "GID"
    #1
    geneinfo <- geneinfo[!duplicated(geneinfo), ]
    geneinfo <- na.omit(geneinfo)
    #2
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
    #3
    if(!is.null(gene2path)){
        if(ncol(gene2path) != 2){
        stop("Dataframe must have only two columns")
        }
        colnames(gene2path) <- c("GID", "PATH")
        gene2path <- gene2path[!duplicated(gene2path), ]
        gene2path <- na.omit(gene2path)
    }else{
        gene2path <- data.frame("GID" = geneinfo$GID,
                        "PATH" = "01100")
    }
    #4
    if(!is.null(gene2symbol)){
      if(ncol(gene2symbol) != 2){
        stop("Dataframe must have only two columns")
      }
      colnames(gene2symbol) <- c("GID", "SYMBOL")
      gene2symbol <- gene2symbol[!duplicated(gene2symbol), ]
      gene2symbol <- na.omit(gene2symbol)
    }else{
      gene2symbol <- data.frame("GID" = geneinfo$GID[1],
                              "SYMBOL" = "")
    }
    #5
    if(!is.null(gene2ensembl)){
      if(ncol(gene2ensembl) != 2){
        stop("Dataframe must have only two columns")
      }
      colnames(gene2ensembl) <- c("GID", "ENSEMBL")
      gene2ensembl <- gene2ensembl[!duplicated(gene2ensembl), ]
      gene2ensembl <- na.omit(gene2ensembl)
    }else{
      gene2ensembl <- data.frame("GID" = geneinfo$GID[1],
                                "ENSEMBL" = "")
    }
    #6
    if(!is.null(gene2refseq)){
      if(ncol(gene2refseq) != 2){
        stop("Dataframe must have only two columns")
      }
      colnames(gene2refseq) <- c("GID", "REFSEQ")
      gene2refseq <- gene2refseq[!duplicated(gene2refseq), ]
      gene2refseq <- na.omit(gene2refseq)
    }else{
      gene2refseq <- data.frame("GID" = geneinfo$GID[1],
                                 "REFSEQ" = "")
    }
    #7
    if(!is.null(gene2pfam)){
      if(ncol(gene2pfam) != 2){
        stop("Dataframe must have only two columns")
      }
      colnames(gene2pfam) <- c("GID", "PFAM")
      gene2pfam <- gene2pfam[!duplicated(gene2pfam), ]
      gene2pfam <- na.omit(gene2pfam)
    }else{
      gene2pfam <- data.frame("GID" = geneinfo$GID[1],
                                "PFAM" = "")
    }
    #8
    if(!is.null(gene2interpro)){
      if(ncol(gene2interpro) != 2){
        stop("Dataframe must have only two columns")
      }
      colnames(gene2interpro) <- c("GID", "INTERPRO")
      gene2interpro <- gene2interpro[!duplicated(gene2interpro), ]
      gene2interpro <- na.omit(gene2interpro)
    }else{
      gene2interpro <- data.frame("GID" = geneinfo$GID[1],
                              "INTERPRO" = "")
    }
    #9
    if(!is.null(gene2reactome)){
      if(ncol(gene2reactome) != 2){
        stop("Dataframe must have only two columns")
      }
      colnames(gene2reactome) <- c("GID", "REACTOME")
      gene2reactome <- gene2reactome[!duplicated(gene2reactome), ]
      gene2reactome <- na.omit(gene2reactome)
    }else{
      gene2reactome <- data.frame("GID" = geneinfo$GID[1],
                                  "REACTOME" = "")
    }
    #10
    if(!is.null(gene2ko)){
      if(ncol(gene2ko) != 2){
        stop("Dataframe must have only two columns")
      }
      colnames(gene2ko) <- c("GID", "KO")
      gene2ko <- gene2ko[!duplicated(gene2ko), ]
      gene2ko <- na.omit(gene2ko)
    }else{
      gene2ko <- data.frame("GID" = geneinfo$GID[1],
                                  "KO" = "")
    }
    #11
    if(!is.null(gene2entrezid)){
      if(ncol(gene2entrezid) != 2){
        stop("Dataframe must have only two columns")
      }
      colnames(gene2entrezid) <- c("GID", "ENTREZID")
      gene2entrezid <- gene2entrezid[!duplicated(gene2entrezid), ]
      gene2entrezid <- na.omit(gene2entrezid)
    }else{
      gene2entrezid <- data.frame("GID" = geneinfo$GID[1],
                            "ENTREZID" = "")
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
    package <- suppressWarnings(makeOrgPackage(
    gene_info = geneinfo,
    symbol = gene2symbol,
    refseq = gene2refseq,
    entrezid = gene2entrezid,
    go = gene2go,
    path = gene2path,
    ko = gene2ko,
    pfam = gene2pfam,
    interpro = gene2interpro,
    reactome = gene2reactome,
    ensembl = gene2ensembl,
    version = version,
    maintainer = maintainer,
    author = author,
    outputDir = outputDir,
    tax_id = tax_id,
    genus = genus,
    species = species,
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
