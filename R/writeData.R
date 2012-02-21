#' Write out the data file for input to sss binary
#'
#' @param data a data frame to be written to file in tempdir()
#' @return the full file path to the file written out
setMethod(
  f = ".writeData",
  signature = "bfrmLinearModel",
  definition = function(object){
    
    tmpMat <- rbind(object@response, object@data)
    fileLoc <- file.path(tempdir(), "data.txt")
    write.table(tmpMat, file=fileLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    object@paramSpec@DataFile <- fileLoc
    object@paramSpec@NContinuousResponses <- 1
    
    return(object)
  }
)

setMethod(
  f = ".writeData",
  signature = "bfrmBinaryModel",
  definition = function(object){
    
    tmpMat <- rbind(object@response, object@data)
    fileLoc <- file.path(tempdir(), "data.txt")
    write.table(tmpMat, file=fileLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    object@paramSpec@DataFile <- fileLoc
    object@paramSpec@NBinaryResponses <- 1
    
    return(object)
  }
)

setMethod(
  f = ".writeData",
  signature = "bfrmCategoricalModel",
  definition = function(object){
    
    tmpMat <- rbind(object@response, object@data)
    fileLoc <- file.path(tempdir(), "data.txt")
    write.table(tmpMat, file=fileLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    object@paramSpec@DataFile <- fileLoc
    object@paramSpec@NCategoricalResponses <- 1
    
    return(object)
  }
)

###########
setMethod(
  f = ".writeData",
  signature = "bfrmSurvivalModel",
  definition = function(object){
    
    tmpMat <- rbind(object@timeToEvent, object@data)
    fileLoc <- file.path(tempdir(), "data.txt")
    write.table(tmpMat, file=fileLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    cenLoc <- file.path(tempdir(), "ResponseMaskFile")
    write.table(object@censor, file=cenLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    object@paramSpec@DataFile <- fileLoc
    object@paramSpec@ResponseMaskFile <- cenLoc
    object@paramSpec@NSurvivalResponses <- 1
    object@paramSpec@NLatentFactors <- 1
    
    return(object)
  }
)


