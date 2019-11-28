# bioAnno
Build Annotation package by using information from KEGG, NCBI, Ensembl
## Description
_bioAnno_ provide wrap functions _fromKEGG_, _fromEnsembl_ and _fromNCBI_ to build annotation package.     
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
##which will build and install package "org.hsa.eg.db" which will include KEGG, GO annotation 
## build Annotation package by using fromEnsembl 
fromEnsembl(species="Homo")
```
## Note
The _bioAnno_ just a package for funs. _bioAnno_ provide wrap function which help me to easily build annotation package.
The _fromAnnotationHub_ function will come soon!!!

## Contact information

For any questions please contact guokai8@gmail.com
