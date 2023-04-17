#' @title build annotation from ensembl
#' @title build annotation from ensembl
#' @importFrom biomaRt useEnsembl listAttributes
#' @importFrom biomaRt useDataset
#' @importFrom biomaRt getBM
#' @importFrom AnnotationForge makeOrgPackage
#' @importFrom utils remove.packages
#' @importFrom utils install.packages
#' @importFrom stats na.omit
#' @param host the ensemble API host,for plant you can use
#'             plants.ensembl.org and for human
#'             and other species you can use uswest.ensembl.org
#' @param species the sepcies you want to search,
#'             you can use listSpecies to get the species name
#' @param anntype the type of function annotation(GO,KEGG,PFAM,InterPro)
#'             you want get from ensemble
#' @param buildall include all prossbile annoation type listed in Ensembl
#' @param author author for the annotation package
#' @param maintainer maintainer for the annotation package
#' @param tax_id taxonomy id for the species
#' @param genus genus name for the annotation package
#' @param version version number for the annotation package
#' @param plant plant or animal species (TRUE/FALSE)
#' @param install install the package or not(default: TRUE)
#' @param pkgname package name you want to choose
#' @param rebuild rebuild the package or not(default: FALSE)
#' @param outputDir temporary output path
#' @examples
#' fromEnsembl(species = "Caenorhabditis elegans", anntype="GO")
#' @author Kai Guo
#' @return annotation package
#' @export
fromEnsembl <- function(species = "Caenorhabditis elegans",
                    host = NULL,
                    anntype = NULL, buildall = TRUE, author = NULL,
                    maintainer = NULL, tax_id = NULL, genus = NULL,
                    version = NULL, plant = FALSE,
                    install = TRUE, pkgname=NULL,outputDir = NULL, rebuild = FALSE){
    if(isTRUE(plant)){
        host = "https://plants.ensembl.org"
        mart = useEnsembl("plants_mart", host = host)
    }else{
     #   mart = useEnsembl("ensembl", mirror = host)
        mart = useEnsembl("ensembl")
        
    }
    dbinfo<-.getmartdb(species, mart)
    if(!is.null(pkgname)){
      dbname1 <- paste0('org.', pkgname, '.eg.db')
    }else{
      dbname1 <- paste0('org.', strsplit(species," ")[[1]][1], '.eg.db')
    }
    if(isTRUE(rebuild)){
        suppressMessages(remove.packages(dbname1))
    }
    if(is_installed(dbname1)){
        suppressMessages(requireNamespace(dbname1, quietly = TRUE))
        cat("You alreay had the annotation package: ", dbname1, " \n")
    }else{
    dbname <- as.character(dbinfo$dbname)
    dataset <- useDataset(dbname, mart = mart)
    attr <- listAttributes(dataset)$name
    if(is.null(anntype)){
        if(isTRUE(buildall)){
        anntype <- c("GO","KEGG","Reactome","PFAM","InterPro")
    }else{
        stop("You need to specific anntation type!\n")
    }
    }
    chr_values <- as.vector(unlist(dbinfo$chr_info$name))
    geneinfo <- getBM(attributes = c("ensembl_gene_id","description"),
                    filters ="chromosome_name", values = chr_values, dataset)
    # geneinfo<-geneinfo[nchar(geneinfo[,2])>1,]
    colnames(geneinfo) <- c("GID", "GENENAME")
    geneinfo <- na.omit(geneinfo)
    gene2symbol <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
            filters = "chromosome_name", values = chr_values, dataset)
    gene2symbol <- gene2symbol[nchar(gene2symbol[,2])>1,]
    colnames(gene2symbol) <- c("GID","SYMBOL")
    gene2entrezid <- getBM(attributes = c("ensembl_gene_id","entrezgene_id"),
            filters ="chromosome_name",values = chr_values, dataset)
    gene2entrezid <- gene2entrezid[nchar(gene2entrezid[,2])>1, ]
    colnames(gene2entrezid) <- c("GID","ENTREZID")
    gene2entrezid <- na.omit(gene2entrezid)
    if(("GO" %in% anntype) & sum(grepl('go_id',attr))>=1){
        gene2go <- getBM(attributes = c("ensembl_gene_id", "go_id",
                    "go_linkage_type"),
                    filters ="chromosome_name", values = chr_values, dataset)
        gene2go <- gene2go[nchar(gene2go[,2])>1, ]
        gene2go <- gene2go[nchar(gene2go[,3])>1, ]
        colnames(gene2go) <- c("GID", "GO", "EVIDENCE")
        gene2go <- na.omit(gene2go)
        gene2go <- gene2go[!duplicated(gene2go), ]
    }else{
        gene2go <- data.frame("GID" = geneinfo$GID,"GO" = "GO:0008150",
        "EVIDENCE" = "IEA")
        cat("Gene Ontology are not list in your annotation database\n")
    }
    if(("KEGG" %in% anntype) & (sum(grepl('kegg_enzyme',attr))>=1)){
        gene2path <- getBM(attributes = c("ensembl_gene_id", "kegg_enzyme"),
                    filters = "chromosome_name", values = chr_values, dataset)
        gene2path[,2] <- sub('\\+.*', '', gene2path[,2])
        gene2path <- gene2path[nchar(gene2path[,2])>1, ]
        colnames(gene2path) <- c("GID", "PATH")
        gene2path <- na.omit(gene2path)
        gene2path <- gene2path[!duplicated(gene2path), ]
    }else{
        gene2path <- data.frame("GID" = geneinfo$GID,"PATH" = "01100")
        cat("KEGG Pathway are not list in your annotation database\n")
    }
    if(("PFAM"%in%anntype)&(sum(grepl('pfam',attr))>=1)){
        gene2pfam <- getBM(attributes = c("ensembl_gene_id","pfam"),
                    filters ="chromosome_name", values = chr_values, dataset)
        gene2pfam <- gene2pfam[nchar(gene2pfam[,2])>1, ]
        colnames(gene2pfam) <- c("GID", "PFAM")
        gene2pfam <- na.omit(gene2pfam)
    }else{
        gene2pfam <- data.frame("GID" = geneinfo$GID, "PFAM" = "PF00001")
        cat("Protein Family are not list in your annotation database\n")
    }
    if(("InterPro" %in% anntype) & (sum(grepl('interpro',attr))>=1)){
        gene2interpro <- getBM(attributes = c("ensembl_gene_id", "interpro"),
            filters = "chromosome_name", values = chr_values, dataset)
        gene2interpro <- gene2interpro[nchar(gene2interpro[, 2])>1, ]
        colnames(gene2interpro) <- c("GID", "INTERPRO")
        gene2interpro <- na.omit(gene2interpro)
    }else{
        gene2interpro <- data.frame("GID" = geneinfo$GID,
            "INTREPRO" = "IPR000001")
        cat("InterPro are not list in your annotation database\n")
    }
    if(("Reactome"%in%anntype) & (sum(grepl('reactome',attr))>=1)){
        if(isTRUE(plant)){
            gene2reactome <- getBM(attributes = c("ensembl_gene_id",
            "plant_reactome_pathway"), filters = "chromosome_name",
            values = chr_values, dataset)
    }else{
        gene2reactome <- getBM(attributes = c("ensembl_gene_id", "reactome"),
        filters = "chromosome_name", values = chr_values, dataset)
    }
        gene2reactome <- gene2reactome[nchar(gene2reactome[,2])>1,]
        colnames(gene2reactome) <- c("GID","REACTOME")
        gene2reactome <- na.omit(gene2reactome)
    }else{
        gene2reactome <- data.frame("GID" = geneinfo$GID,
        "REACTOME" = "RSA0000000")
        cat("Reactome Pathway are not list in your annotation database\n")
    }
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
    if(!is.null(pkgname)){
      species <- pkgname
    }
    package <- suppressWarnings(makeOrgPackage(gene_info = geneinfo,
        symbol = gene2symbol,
        entrezid = gene2entrezid,
        go = gene2go,
        path = gene2path,
        pfam = gene2pfam,
        interpro = gene2interpro,
        reactome = gene2reactome,
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
}
##' @title list species available in Ensembl
##' @importFrom  biomaRt useMart
##' @importFrom  biomaRt listDatasets
##' @param host Ensembl host site
##' @param plant use plant database or not (default: FALSE)
##' @examples
##' listSpecies()
##' @author Kai Guo
##' @return data.frame with species information
##' @export
listSpecies <- function(host = "www", plant = FALSE){
    cat("You could choose different host to get high speed!\n")
    if(isTRUE(plant)){
        host <- "https://plants.ensembl.org"
        mart <- useEnsembl("plants_mart", host = host)
    }else{
        cat("host: 'www', 'uswest', 'useast', 'asia'\n" )
      #  mart <- useEnsembl("ENSEMBL_MART_ENSEMBL", mirror = host)
        mart <- useEnsembl("ensembl")
        
    }
    res <- tryCatch(
        expr = { 
            listDatasets(mart)
        },
        error = function(e){
            NULL
        },
        warning = function(w){
            NULL
        })
    if(is.data.frame(res)) colnames(res)[2] <- "species"
    res
}

