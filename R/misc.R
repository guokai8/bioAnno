##' uppercase the first letter
##' @param x string
##' @return character with first letter uppercase
##' @author Kai Guo
simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep = "", collapse = " ")
}

#' extract GO information from NCBI and filter by taxid
#' @importFrom data.table fread
#' @importFrom R.utils gunzip
#' @importFrom utils download.file
#' @importFrom data.table ":="
#' @param taxid taxonomy id for the species
#' @param species species name(common name,kegg.species.code or scientifc name)
#' @return dataframe with gene2go information
#' @author Kai Guo
.extratGO <- function(taxid = NULL, species = NULL){
    # temp file
    if(is.null(taxid)){
    taxid <- .get.species.info(species)['tax.id']
    }
    tmp<-paste(
    tempfile(),
    "gz",
    sep = "."
    )
    # import Gene to Gene Ontology from NCBI Gene database
    download.file(
    "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz",
    quiet = TRUE,
    destfile = tmp
    )
    # uncompress
    gunzip(tmp)
    # read the file (linux and windows)
    gene2go = fread(
    sub("\\.gz", "", tmp),
    verbose = FALSE,
    showProgress = FALSE
    )
    # select columns and rename
    gene2go <- unique(
    gene2go[,c(seq_len(4)), with = FALSE]
    )
    colnames(gene2go) <- c("taxid", "GID", "GO", "EVIDENCE")

    # convert columns in character
    gene2go[,
        `:=`(
        taxid = as.character(gene2go$taxid),
        GID = as.character(gene2go$GID)
        )
        ]
    # filter with taxid
    gene2go <- as.data.frame(gene2go)
    gene2go <- gene2go[gene2go$taxid == taxid,2:4]
    return(gene2go)
}


#' @importFrom data.table fread
#' @importFrom R.utils gunzip
#' @importFrom utils download.file
#' @importFrom data.table ":="
#' @param taxid taxonomy id for the species
#' @param species species name(common name,kegg.species.code or scientifc name)
#' @author Kai Guo
.extratGene <- function(taxid = NULL, species = NULL){
    if(is.null(taxid)){
    taxid <- .get.species.info(species)['tax.id']
    }
    # temp file
    tmp <- paste(
    tempfile(),
    "gz",
    sep="."
    )
    # import Gene to Gene Ontology from NCBI Gene database
    download.file(
    "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz",
    quiet = TRUE,
    destfile = tmp
    )
    # uncompress
    gunzip(tmp)
    # read the file (linux and windows)
    gene2info <- fread(
    sub("\\.gz", "", tmp),
    verbose=FALSE,
    showProgress=FALSE
    )
    # select columns and rename
    gene2info <- unique(
    gene2info[, c(seq_len(3), 9), with = FALSE]
    )
    colnames(gene2info) <- c("taxid", "GID", "SYMBOL", "GENENAME")

    # convert columns in character
    gene2info[,
        `:=`(
            taxid = as.character(gene2go$taxid),
            GID = as.character(gene2go$GID)
        )
        ]
    # filter with taxid
    gene2go <- gene2go[taxid == taxid,2:4]
    return(as.data.frame(gene2info))
}
##' modified from pathview kegg.species.code
##' @importFrom utils data
#' @param species species name(common name,kegg.species.code or scientifc name)
#' @param na.rm TRUE/FALSE
#' @return character with species name
#' @author Kai Guo
.get.species.info <- function (species = "hsa", na.rm = FALSE)
{
    nspec <- length(species)
    if(!exists("korg")) data(korg)
    ridx <- match(species, korg[, seq_len(5)]) %% nrow(korg)
    nai <- is.na(ridx)
    if (sum(nai) > 0) {
    na.msg <- sprintf("Unknown species '%s'!", paste(species[nai],
        sep = "", collapse = "', '"))
        message("Note: ", na.msg)
    }
    if (sum(nai) == nspec) {
        stop.msg <- "All species are invalid!"
        stop(stop.msg)
    }
    if (any(ridx[!nai] == 0)) ridx[!nai & ridx == 0] <- nrow(korg)
    if (isTRUE(na.rm)) ridx = ridx[!nai]
    species.info <- korg[ridx, ]
    return(species.info)
}

##' @title get species information in Ensembl
##' @importFrom dplyr select_
##' @importFrom dplyr collect
##' @importFrom dplyr pull
##' @importFrom magrittr %>%
##' @importFrom biomaRt listDatasets
##' @importFrom jsonlite fromJSON
##' @param species species
##' @param mart biomaRt mart
##' @return list with species information
##' @author Kai Guo
.getmartdb <- function(species, mart){
    lhs <- listDatasets(mart)
    spe <- simpleCap(species);
    spe <- gsub(' ', '\\\\s', spe)
    sel <- grepl(spe, lhs$description, ignore.case = TRUE)
    tmp <- lhs[sel,]
    dataset <- tmp%>%select_(~dataset)%>%collect%>%pull(1)
    if((length(dataset) == 0) | (length(dataset) > 1)){
        stop("Maybe you need first check the avaliable database by
        using listSpecies()\n")
    }
    chr <- tmp%>%select_(~description)%>%collect%>%pull(1)
    organism <- gsub(' ', '_', sub(' genes.*', '', chr))
    if(organism == "Oryza_sativa_Japonica"){
        organism = "Oryza_sativa"
    }
    if(mart@biomart =="plants_mart"){
    pre_site <- "http://rest.ensembl.org/info/assembly/"
    }else{
    pre_site <- "http://rest.ensembl.org/info/assembly/"
    }
    tryCatch({
    chr_d <-
        fromJSON(
        paste0(
        pre_site,
        organism,
        "?content-type=application/json"
        )
    )
    }, error = function(e)
    stop(
    "The API 'http://rest.ensembl.org' does not seem to work properly.
    Are you connected to the internet? Is the homepage
    'http://rest.ensembl.org' currently available?", call. = FALSE
    ))
    chr_info <- chr_d$top_level_region
    chr_version <- chr_d$assembly_name
    chr_assembly_date <- chr_d$assembly_date
    rhs <- list(dbname = dataset, chr_info = chr_info,
        chr_version = chr_version, chr_assembly_date = chr_assembly_date)
    return(rhs)
}

#'@title get database name by using species name
#'@param species species name
#'@return character with database name
#' @author Kai Guo
.getdbname <- function(species = species){
    species = tryCatch(match.arg(species,c("anopheles",
        "arabidopsis", "bovine", "celegans", "canine", "fly", "zebrafish",
        "ecoli", "ecsakai", "chicken", "human", "mouse", "rhesus", "malaria",
        "chipm", "rat",
        "toxoplasma", "streptomyces", "pig", "yeast", "xenopus", "warm")),
        error=function(cond){return("unsupported")})
    if (species == "anopheles") {
    dbname <- "org.Ag.eg.db.sqlite"
    } else if (species == "bovine") {
    dbname <- "org.Bt.eg.db.sqlite"
    } else if (species == "canine") {
    dbname <- "org.Cf.eg.db.sqlite"
    } else if (species == "worm" || species == "celegans") {
    dbname <- "org.Ce.eg.db.sqlite"
    } else if (species == "chicken") {
    dbname <- "org.Gg.eg.db.sqlite"
    } else if (species == "ecolik12") {
    dbname <- "org.EcK12.eg.db.sqlite"
    } else if (species == "ecsakai") {
    dbname <- "org.EcSakai.eg.db.sqlite"
    } else if (species == "fly") {
    dbname <- "org.Dm.eg.db.sqlite"
    } else if (species == "human") {
    dbname <- "org.Hs.eg.db.sqlite"
    } else if (species == "chipm") {
    dbname <- "org.Pt.eg.db.sqlite"
    }else if (species == "mouse") {
    dbname <- "org.Mm.eg.db.sqlite"
    } else if (species == "pig") {
    dbname <- "org.Ss.eg.db.sqlite"
    } else if (species == "rat") {
    dbname <- "org.Rn.eg.db.sqlite"
    } else if (species == "rhesus") {
    dbname <- "org.Mmu.eg.db.sqlite"
    } else if (species == "xenopus") {
    dbname <- "org.Xl.eg.db.sqlite"
    } else if (species == "zebrafish") {
    dbname <- "org.Dr.eg.db.sqlite"
    } else {
    dbname <- NULL
    }
    return(dbname)
}
##' @author Kai Guo
readidx <- function()
{
    n <- readline(prompt = "Enter an index: ")
    return(as.integer(n))
}
##' @title check package installed or not
##' @param pkg package name
##' @return TRUE/FALSE
##' @author Kai Guo
is_installed <- function(pkg) {
    nzchar(system.file(package = pkg))
}
##' @title show the package path
##' @param package the full path of the package
##' @return whole path for the package
##' @author Kai Guo
.show.path <- function(package){
    cat("################################################################\n")
    cat("Please find your annotation package in ...\n")
    cat(package,"\n")
    cat("You can install it by using\n")
    cat(paste0("install.packages(\"",package,",repos = NULL,type='source')"),"\n")
    cat("################################################################\n")
}
##' @title show the package content
##' @importFrom RSQLite dbConnect
##' @importFrom RSQLite SQLite
##' @importFrom RSQLite dbListTables
##' @importFrom RSQLite dbReadTable
##' @importFrom RSQLite dbDisconnect
##' @param package the full path of the package
##' @return vector 
##' @author Kai Guo
.show.tables <- function(package){
    pkg <- basename(package)
    path <- paste0(package,"/inst/extdata/",sub('.db','.sqlite',pkg))
    con <- dbConnect(SQLite(),path)
    cat("Here are the tables in the package",pkg,"...\n")
    dblist<-dbListTables(con)
    cat(dblist,"\n")
    cat("################################################################\n")
    dbDisconnect(con)
}
##' @title get annotataion table from temporary package
##' @importFrom RSQLite dbConnect
##' @importFrom RSQLite SQLite
##' @importFrom RSQLite dbReadTable
##' @importFrom RSQLite dbDisconnect
##' @importFrom dplyr left_join
##' @param path full path for the temporary package
##' @param table a character  indicate the table you want extract
##' @examples
##' data(ath)
##' pack <- fromOwn(geneinfo = ath, install = FALSE, species ="test")
##' # head(getTable(path = pack, table = "gene_info"))
##' @export
##' @return data.frame 
##' @author Kai Guo

getTable <- function(path,table="go_all"){
    pkg <- basename(path = path)
    path <- paste0(path,"/inst/extdata/",sub('.db','.sqlite',pkg))
    con <- dbConnect(SQLite(),path)
    gene_info <- dbReadTable(con,"gene_info")
    anno <- dbReadTable(con,table)
    res <- left_join(gene_info, anno, by=c("X_id"="X_id"))
    colnames(res)[1] <- "ID"
    dbDisconnect(con)
    return(res)
}

