# bioAnno <a href="https://travis-ci.org/guokai8/bioAnno"><img src="https://travis-ci.org/guokai8/bioAnno.svg" alt="Build status"></a>    [![](https://img.shields.io/badge/devel%20version-0.0.6-green.svg)](https://github.com/guokai8/bioAnno)

Build Annotation package by using information from __KEGG__, __NCBI__, __Ensembl__ and return OrgDb object such as org.Hs.eg.db. The _bioAnno_ package support all organisms list in __Ensembl__, __KEGG__, __NCBI__.  
## Description
_bioAnno_ provide wrap functions _fromKEGG_, _fromEnsembl_,_fromNCBI_ and _fromAnnoHub_ to build annotation package. 
KEGG species code is suggested to use except _fromEnsembl_ which require scientific name.
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
## Note
The _bioAnno_ just a package for funs. _bioAnno_ provide wrap function which help me to easily build annotation package.
The _fromAnnoHub_ function is ready to use. 

## Contact information

For any questions please contact guokai8@gmail.com
