is.MRFA <- function(object){
  if (inherits(object, "MRFA")) return(TRUE)
  else return(FALSE)
}
