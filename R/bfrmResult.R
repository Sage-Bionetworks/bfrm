## METHOD TO READ IN SUMMARY OUTPUT FROM bfrm FOR DIFFERENT MODEL TYPES
#####


setMethod(
  f = ".readResult",
  signature = c("bfrmModel", "character"),
  definition = function(object, outputDir){
    
    theseFiles <- as.list(sub(".txt", "", list.files(outputDir, pattern=".txt"), fixed=T))
    
    outList <- lapply(theseFiles, function(x){
      obj <- as.matrix(read.delim(file.path(outputDir, paste(x, ".txt", sep="")), header=F, as.is=T, strip.white=T))
      tmp <- colSums(is.na(obj))
      obj <- obj[, !(tmp==dim(obj)[1])]
      
#       if(x=="mA"){
# #        rownames(obj) ## y and x in that order
# #        colnames(obj)
#       } else if(x=="mF"){
# #        rownames(obj)
#         if(!is.null(colnames(object@data)))
#           colnames(obj) <- colnames(object@data)
#       } else if(x=="mPib"){
# #        rownames(obj)
# #        colnames(obj)
#       } else if(x=="mPostPib"){
# #        rownames(obj) ## y and x in that order
# #        colnames(obj) ## SAME AS mA
#       } else if(x=="m"){
# #        rownames(obj)
# #        colnames(obj) ## y and x in that order
#       } else if(x=="mTau"){
# #        rownames(obj)
# #        colnames(obj) ## SAME AS mA
#       } else if(x=="mPsi"){
# #        rownames(obj)
# #        colnames(obj) ## y and x in that order
#       } else if(x=="mVariablesIn"){ ## ONLY IN EVOLUTIONARY MODE
# #        rownames(obj)
# #        colnames(obj)
#       } else if(x=="mX"){
#         if(!is.null(rownames(object@data)))
#           rownames(obj) <- rownames(object@data)
#         if(!is.null(colnames(object@data)))
#           colnames(obj) <- colnames(object@data)
#       } else if(x=="mZ"){
# #        rownames(obj)
#         if(!is.null(colnames(object@y)))
#           colnames(obj) <- colnames(object@y)
#       }
      
      return(obj)
    })
    names(outList) <- theseFiles
    
    thisClass <- sub("Model", "Result", class(object))
    
    newObj <- new(thisClass,
                  model = object,
                  results = outList)
    
    return(newObj)
  }
)

