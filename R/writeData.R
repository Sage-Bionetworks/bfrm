
setMethod(
  f = ".writeData",
  signature = "bfrmModel",
  definition = function(object){
    
    tmpMat <- rbind(object@y, object@data)
    fileLoc <- tempfile(pattern="data", tmpdir=tempdir(), fileext=".txt")
    write.table(tmpMat, file=fileLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    tmpMask <- ifelse(is.na(object@data), 1, 0)
    maskLoc <- tempfile(pattern="xmaskfile", tmpdir=tempdir(), fileext=".txt")
    write.table(tmpMask, file=maskLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    object@paramSpec@xmaskfile <- maskLoc
    object@paramSpec@datafile <- fileLoc
    object@paramSpec@nobservations <- ncol(tmpMat)
    object@paramSpec@nvariables <- nrow(tmpMat)
    
    if( (length(object@design) + length(object@control)) > 0 ){
      tmpH <- rbind(1, object@design, object@control)
      tmpH <- t(tmpH)
      hLoc <- tempfile(pattern="hfile", tmpdir=tempdir(), fileext=".txt")
      write.table(tmpH, file=hLoc, sep="\t", quote=F, row.names=F, col.names=F)
      object@paramSpec@hfile <- hLoc
      object@paramSpec@ndesignvariables <- nrow(object@design) + 1
      object@paramSpec@ncontrolvariables <- nrow(object@control)
    }
    
    if( length(object@ymask) > 0 ){
      cenLoc <- tempfile(pattern="responseMaskFile", tmpdir=tempdir(), fileext=".txt")
      write.table(object@ymask, file=cenLoc, sep="\t", quote=F, row.names=F, col.names=F)
      object@paramSpec@responsemaskfile <- cenLoc
    } else{
      object@paramSpec@responsemaskfile <- "NA"
    }
    
    return(object)
  }
)

