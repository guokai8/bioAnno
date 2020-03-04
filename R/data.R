##' @title korg
##' @name korg
##' @description korg include species information from KEGG database.
##'              korg data was modified from
##'              (https://pathview.uncc.edu/data/korg.tsv)
##' @format A matrix with five columns:
##' \describe{
##' \item{ktax.id}{KEGG taxonomy ID}
##' \item{tax.id}{NCBI taxonomy ID}
##' \item{kegg.code}{KEGG species code}
##' \item{scientific.name}{Scientific name of species}
##' \item{common.name}{common name of species}
##' }
##' @examples
##' \dontrun{
##' korg
##' }
##'
"korg"

##' @title TAIR10 geneid to ENTREZID
##' @name ath
##' @description ath include TAIR10 ID with matched ENTREZID information
##' @format A matrix with two columns:
##' \describe{
##' \item{GID}{TAIR GENE ID}
##' \item{ENTREZID}{NCBI ENTREZID ID}
##' }
##' @examples
##' \dontrun{
##' ath
##' }
##'
"ath"
