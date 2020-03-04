# bioAnno <a href="https://travis-ci.org/guokai8/bioAnno"><img src="https://travis-ci.org/guokai8/bioAnno.svg" alt="Build status"></a>  [![Project Status:](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![](https://img.shields.io/badge/devel%20version-0.0.6-green.svg)](https://github.com/guokai8/bioAnno)    

Build Annotation package by using information from __KEGG__, __NCBI__, __Ensembl__ and return OrgDb object such as org.Hs.eg.db. The _bioAnno_ package support all organisms list in __Ensembl__, __KEGG__, __NCBI__.  
## Description
With the increasing of high throughput data generated, the requirement for
having annotation package ready is necessary for people doing functional 
enrichment analysis, id conversion and other type related analysis.
_bioAnno_ provide wrap functions include _fromKEGG_, _fromEnsembl_, 
_fromNCBI_ and _fromAnnoHub_ to build annotation package. 
And you can easily to build annotation package with 
the KEGG species code (except _fromEnsembl_ which require scientific name).
## Installation
```
library(devtools)
install_github("guokai8/bioAnno")
``` 

## Software Usage

```
library(bioAnno)
## build Annotation package by using fromKEGG
fromKEGG(species="hsa")
## which will build and install package "org.hsa.eg.db" which will include KEGG, GO annotation 
## build Annotation package by using fromEnsembl 
fromEnsembl(species="Human")
## build from AnnotationHub
fromAnnHub(species="human")
```
### 2.4 Main Functions
--  _fromKEGG_ build annotation package by extracting annotation information 
    from Kyoto Encyclopedia of Genes and Genomes database (KEGG). 
    You can use kegg species code as query species name.

-- _fromNCBI_ build annotation package by extracting annotation information from
    NCBI database.

-- _fromENSEMBL_ build annotation package by extracting annotation information 
    fromENSEMBL database. It includes function to build annotaion package for 
    plant with parameter plant = TRUE
-- _fromAnnhub_ build annotation package with the AnnotationHub package 
## Note
The _bioAnno_ provide wrap function which help me to easily build annotation package.
The _fromAnnoHub_ function is ready to use. 

## Contact information

For any questions please contact guokai8@gmail.com
