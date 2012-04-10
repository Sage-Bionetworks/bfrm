
setMethod(
  f = ".writeData",
  signature = "bfrmModel",
  definition = function(object){
    
    fileLoc <- tempfile(pattern="data", tmpdir=tempdir(), fileext=".txt")
    write.table(object@data, file=fileLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    tmpMask <- ifelse(is.na(object@data), 1, 0)
    maskLoc <- tempfile(pattern="xmaskfile", tmpdir=tempdir(), fileext=".txt")
    write.table(tmpMask, file=maskLoc, sep="\t", quote=F, row.names=F, col.names=F)
    
    object@paramSpec@xmaskfile <- maskLoc
    object@paramSpec@datafile <- fileLoc
    object@paramSpec@nobservations <- ncol(object@data)
    object@paramSpec@nvariables <- nrow(object@data)
    
    if( (length(object@design) + length(object@control)) > 0 ){
      tmpH <- rbind(1, object@design, object@control)
      tmpH <- t(tmpH)
      hLoc <- tempfile(pattern="hfile", tmpdir=tempdir(), fileext=".txt")
      write.table(tmpH, file=hLoc, sep="\t", quote=F, row.names=F, col.names=F)
      object@paramSpec@hfile <- hLoc
      object@paramSpec@ndesignvariables <- nrow(object@design) + 1
      object@paramSpec@ncontrolvariables <- nrow(object@control)
    }
    
    return(object)
  }
)

