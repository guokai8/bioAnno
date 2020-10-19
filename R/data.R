##' @title korg
##' @name korg
##' @description korg include species information from KEGG database.
##'              korg data was modified from
##'              (https://pathview.uncc.edu/data/korg.tsv)
##' @format A matrix with five columns:
##' \describe{
##' \item{ktax.id}{the KEGG taxonomy ID}
##' \item{tax.id}{the NCBI taxonomy ID}
##' \item{kegg.code}{the KEGG species code}
##' \item{scientific.name}{Scientific name of species}
##' \item{common.name}{common name of species}
##' }
##' @examples
##' head(korg)
##'
"korg"

##' @title TAIR10 geneid to ENTREZID
##' @name ath
##' @description The 'ath' dataset include the annotation information collected
##' form the TAIR10 database(htps://arabidopsis.org/download/index-auto.jsp
##' %3Fdir%3D%252Fdownload_files%252FGenes%252FTAIR10_genome_release). 
##' @format A data.frame with two columns: 
##' \describe{
##' \item{GID}{The arabidopsis GENE ID}
##' \item{ENTREZID}{NCBI ENTREZID ID for the arabidopsis}
##' }
##' @examples
##' head(ath)
##'
"ath"
