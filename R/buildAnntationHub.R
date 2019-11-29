#' extract annotation database by using AnnotationHub
#' @importFrom AnnotationHub query
#' @author Kai Guo
#' @export
fromAnnHub<-function(species,useEnsembl=FALSE,useNCBI=FALSE,useUCSC=FALSE){
  ah<-AnnotationHub()
  if(isTRUE(useEnsembl)){
    ah <- query(ah, "EnsDb")
  }
  else if(isTRUE(useNCBI)){
    ah <- query(ah,pattern = c(species,"NCBI"))
  }else{
    ah <- query(ah,pattern = species)
  }
  #idx <- grep('*eg')
}




