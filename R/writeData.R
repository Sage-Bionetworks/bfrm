#' Write out the data file for input to sss binary
#'
#' @param data a data frame to be written to file in tempdir()
#' @return the full file path to the file written out
setMethod(
  f = ".writeData",
  signature = "bfrmLinearModel",
  definition = function(object){
    
    tmpMat <- rbind(object@response, object@data)
    fileLoc <- tempfile(pattern="data", tmpdir=tempdir(), fileext=".txt")
    write.table(tmpMat, file=fileLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    tmpMask <- is.na(object@data)
    maskLoc <- tempfile(pattern="xmaskfile", tmpdir=tempdir(), fileext=".txt")
    write.table(tmpMask, file=maskLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    object@paramSpec@xmaskfile <- maskLoc
    object@paramSpec@datafile <- fileLoc
    object@paramSpec@ncontinuousresponses <- 1
    object@paramSpec@nobservations <- ncol(tmpMat)
    object@paramSpec@nvariables <- nrow(tmpMat)
    
    return(object)
  }
)

setMethod(
  f = ".writeData",
  signature = "bfrmBinaryModel",
  definition = function(object){
    
    tmpMat <- rbind(object@response, object@data)
    fileLoc <- tempfile(pattern="data", tmpdir=tempdir(), fileext=".txt")
    write.table(tmpMat, file=fileLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    tmpMask <- is.na(object@data)
    maskLoc <- tempfile(pattern="xmaskfile", tmpdir=tempdir(), fileext=".txt")
    write.table(tmpMask, file=maskLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    object@paramSpec@xmaskfile <- maskLoc
    object@paramSpec@datafile <- fileLoc
    object@paramSpec@nbinaryresponses <- 1
    object@paramSpec@nobservations <- ncol(tmpMat)
    object@paramSpec@nvariables <- nrow(tmpMat)
    
    return(object)
  }
)

setMethod(
  f = ".writeData",
  signature = "bfrmCategoricalModel",
  definition = function(object){
    
    tmpMat <- rbind(object@response, object@data)
    fileLoc <- tempfile(pattern="data", tmpdir=tempdir(), fileext=".txt")
    write.table(tmpMat, file=fileLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    tmpMask <- is.na(object@data)
    maskLoc <- tempfile(pattern="xmaskfile", tmpdir=tempdir(), fileext=".txt")
    write.table(tmpMask, file=maskLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    object@paramSpec@xmaskfile <- maskLoc
    object@paramSpec@DataFile <- fileLoc
    object@paramSpec@ncategoricalresponses <- 1
    object@paramSpec@nobservations <- ncol(tmpMat)
    object@paramSpec@nvariables <- nrow(tmpMat)
    
    return(object)
  }
)

###########
setMethod(
  f = ".writeData",
  signature = "bfrmSurvivalModel",
  definition = function(object){
    
    tmpMat <- rbind(object@timeToEvent, object@data)
    fileLoc <- tempfile(pattern="data", tmpdir=tempdir(), fileext=".txt")
    write.table(tmpMat, file=fileLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    ## HACK IN GOOFY NOTATION FOR BFRM CODE
    ## 0 IS AN EVENT - AND 2 IS A CENSOR
    tmpCen <- ifelse(object@censor==1, 0, 2)
    cenLoc <- tempfile(pattern="responseMaskFile", tmpdir=tempdir(), fileext=".txt")
    write.table(tmpCen, file=cenLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    tmpMask <- is.na(object@data)
    maskLoc <- tempfile(pattern="xmaskfile", tmpdir=tempdir(), fileext=".txt")
    write.table(tmpMask, file=maskLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    object@paramSpec@xmaskfile <- maskLoc
    object@paramSpec@datafile <- fileLoc
    object@paramSpec@responsemaskfile <- cenLoc
    object@paramSpec@nsurvivalresponses <- 1
    object@paramSpec@nlatentfactors <- 1
    object@paramSpec@nobservations <- ncol(tmpMat)
    object@paramSpec@nvariables <- nrow(tmpMat)
    
    return(object)
  }
)


