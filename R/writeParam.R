#' Write out the param file for input to sss binary
#'
#' @param paramSpec a list to be written to file in tempdir()
#' @return the full file path to the file written out
setMethod(
  f = ".writeParam",
  signature = "bfrmParam",
  definition = function(paramSpec){
    
    thisOut <- sapply(slotNames(paramSpec), function(x){
      paste(x, " = ", slot(paramSpec, x), sep="")
    })
                      
    fileLoc <- file.path(tempdir(), "parameter.txt")
    write.table(thisOut, file=fileLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    return(fileLoc)
  }
)



