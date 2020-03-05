#' @title extract annotation database by using AnnotationHub
#' @importFrom AnnotationHub AnnotationHub
#' @importFrom AnnotationHub query
#' @importFrom AnnotationDbi select
#' @importFrom AnnotationDbi keys
#' @importFrom AnnotationDbi columns
#' @importFrom AnnotationForge makeOrgPackage
#' @importFrom RSQLite dbGetQuery
#' @importFrom utils remove.packages
#' @importFrom utils install.packages
#' @importFrom stats na.omit
#' @param species species name(common name,kegg.species.code or scientifc name)
#' @param author author for the annotation package
#' @param maintainer maintainer for the annotation package
#' @param tax_id taxonomy id for the species
#' @param genus genus for the species
#' @param version version for the annotation package
#' @param install install the package or not(default: TRUE)
#' @param outputDir temporary file path
#' @param rebuild rebuild the package or not(default: FALSE)
#' @examples
#' \donttest{
#' fromAnnHub(species="ath")
#' }
#' @author Kai Guo
#' @return annotation package
#' @export
fromAnnHub<-function(species, author = NULL,
        maintainer = NULL, tax_id = NULL, genus = NULL,
        version = NULL,
        install = TRUE, outputDir = NULL, rebuild = FALSE){
    dbi <- .getdbname(species)

    if(is.null(dbi)){
        dbi <- .get.species.info(species = species)
        dbi <- dbi["scientific.name"]
        dbi <- paste0(unlist(strsplit(dbi, ' '))[c(1,2)], collapse = " ")
    }
    species <- gsub(' .*', '', species)
    dbname <- paste0("org.", species, ".eg.db")
    if(isTRUE(rebuild)){
        suppressMessages(remove.packages(dbname))
    }
    if(is_installed(dbname)){
        suppressMessages(requireNamespace(dbname, quietly = TRUE))
        cat("You alreay had the annotation package: ", dbname, " \n")
    }else{
    ah <- AnnotationHub()
    ah <- query(ah, dbi)
    ahdb <- ah$title
    names(ahdb) <- ah$ah_id
    ### only need ors.xxx.eg.xxx or org.xxx.db.sqlite
    idx<-grep('^org', grep('\\.[eg|db]', ahdb,value = TRUE), value = TRUE)
    #idx <- grep(sub(' .*','',sub(' ','_',dbi)),ahdb,value=T)
    if(length(idx)>1){
        cat("Please select which database you want to use (1,2,3,...): \n")
        for(i in seq_len(length(idx))){
        cat(i, ":", idx[i], "\n")
    }
    idd <- readidx()
    idx <- idx[idd]
    }
    idn <- names(idx)
    res <- ah[[idn]]
    packinfo <- dbGetQuery(res$conn, "select * from metadata;")
    geneinfo <- select(res,keys = keys(res), columns=c("GENENAME"))
    geneinfo <- na.omit(geneinfo)
    colnames(geneinfo)[1] <- "GID"
    gene2refseq <- select(res, keys = keys(res), columns = c("REFSEQ"))
    gene2refseq <- na.omit(gene2refseq)
    colnames(gene2refseq)[1] <- "GID"
    gene2symbol <- select(res, keys = keys(res), columns = c("SYMBOL"))
    gene2symbol <- na.omit(gene2symbol)
    colnames(gene2symbol)[1] <- "GID"
    gene2go <- select(res, keys = keys(res), columns = c("GOALL","EVIDENCEALL"))
    gene2go <- gene2go[, c(1,2,3)]
    gene2go <- gene2go[!duplicated(gene2go),]
    gene2go <- na.omit(gene2go)
    colnames(gene2go) <- c("GID", "GO", "EVIDENCE")
    pathway = FALSE
    if("PATH" %in% columns(res)){
        gene2path <- select(res, keys = keys(res), columns = c("PATH"))
        colnames(gene2path) <- c("GID", "PATH")
        gene2path <- gene2path[!duplicated(gene2path), ]
        gene2path <- na.omit(gene2path)
        pathway <- TRUE
    }
    if(is.null(version)){
        version <- "0.0.1"
    }
    if(is.null(tax_id)){
        tax_id <- packinfo[6,2]
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
        species <- species
    }
    if(is.null(outputDir)){
        outputDir <- tempdir()
    }
    if(isTRUE(install)){
        if(isTRUE(pathway)){
        package<-makeOrgPackage(
            gene_info = geneinfo,
            refseq = gene2refseq,
            symbol = gene2symbol,
            go = gene2go,
            path = gene2path,
            maintainer = maintainer,
            author = author,
            outputDir = outputDir,
            tax_id = "tax_id",
            genus = "",
            species = species,
            version = version,
            verbose = FALSE,
            goTable = "go")
    }else{
        package<-makeOrgPackage(
        gene_info = geneinfo,
        refseq = gene2refseq,
        symbol = gene2symbol,
        go = gene2go,
        maintainer = maintainer,
        author = author,
        outputDir = outputDir,
        tax_id = "tax_id",
        genus = "",
        species = species,
        version = version,
        verbose = FALSE,
        goTable = "go")
    }
    install.packages(package, repos = NULL, type = "source")
    unlink(package, recursive = TRUE)
    }else{
        return(res)
    }
    }
}
