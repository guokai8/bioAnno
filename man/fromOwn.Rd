% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildOwn.R
\name{fromOwn}
\alias{fromOwn}
\title{build Own annotation database with user defined annotation file}
\usage{
fromOwn(
  geneinfo = geneinfo,
  keytype = NULL,
  gene2go = NULL,
  gene2path = NULL,
  gene2symbol = NULL,
  gene2refseq = NULL,
  gene2ensembl = NULL,
  gene2pfam = NULL,
  gene2reactome = NULL,
  gene2ko = NULL,
  gene2interpro = NULL,
  gene2entrezid = NULL,
  gene2biocyc = NULL,
  gene2kd = NULL,
  gene2fundo = NULL,
  gene2gad = NULL,
  version = NULL,
  maintainer = NULL,
  author = NULL,
  outputDir = NULL,
  tax_id = NULL,
  genus = NULL,
  species = NULL,
  install = TRUE,
  pkgname = NULL,
  rebuild = FALSE
)
}
\arguments{
\item{geneinfo}{gene information table with two columns
as default("GID","DESCRIPTION")}

\item{keytype}{key type for building the annotation db}

\item{gene2go}{Gene Onotoly information for  genes}

\item{gene2path}{KEGG Pathway information for genes}

\item{gene2symbol}{SYMBOL information for genes}

\item{gene2refseq}{REFSEQ or KO information for genes}

\item{gene2ensembl}{ENSEMBL or KO information for genes}

\item{gene2pfam}{PFAM information for genes}

\item{gene2reactome}{REACTOME Pathway or KO information for genes}

\item{gene2ko}{KO information for genes}

\item{gene2interpro}{INTERPRO information for genes}

\item{gene2entrezid}{ENTREZID information for genes}

\item{gene2biocyc}{BIOCYC information for genes}

\item{gene2kd}{KEGG DISEASE information for genes}

\item{gene2fundo}{FunDO information for genes}

\item{gene2gad}{GAD information for genes}

\item{version}{version for the annotation package}

\item{maintainer}{maintainer for the annotation package}

\item{author}{author for the annotation package}

\item{outputDir}{temporary output path}

\item{tax_id}{taxonomy id for the species}

\item{genus}{genus for the species}

\item{species}{species name(common name,kegg.species.code or scientifc name)}

\item{install}{install the package or not(default: TRUE)}

\item{pkgname}{package name you want to choose}
}
\value{
annotation package
}
\description{
build Own annotation database with user defined annotation file
}
\examples{
## build your own annotation for Arabidopsis thaliana
data(ath)
fromOwn(geneinfo = ath, install = FALSE)
}
\author{
Kai Guo
}
