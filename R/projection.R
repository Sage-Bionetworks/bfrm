## PROJECT FACTORS ONTO ANOTHER DATASET
## SUMMARIZATION OF RESULTS AS PER MATLAB CODE SENT FROM JOE LUCAS
#####

setMethod(
  f = "projection",
  signature = c("bfrmResult", "matrix"),
  definition = function(object, newdata){
    
    ## CHECK IF EITHER SETS OF ROWNAMES ARE NULL - IF SO, MUST ASSUME ORDERED SAME
    odr <- is.null(rownames(object@model@data))
    ndr <- is.null(rownames(newdata))
    if(odr | ndr){
      warning("rows not named in newdata or original data - must be in same order for projections to be accurate")
      pidNew <- 1:nrow(newdata)
      pidOrig <- 1:nrow(object@model@data)
    } else{
      pidNew <- rownames(newdata)
      pidOrig <- rownames(object@model@data)
    }
    
    ## SUBSET VARIABLES IF NECESSARY
    if( any(names(object@results) == "mVariablesIn") ){
      pidOrig <- pidOrig[as.numeric(object@results$mVariablesIn)]
      pidNew <- pidNew[as.numeric(object@results$mVariablesIn)]
    }
    
    if( ! all(pidOrig %in% pidNew) ){
      stop("newdata must contain all features present in original data")
    }
    
    useMat <- newdata[pidOrig, ]
    
    B <- object@results$mA * object@results$mPostPib
    mpi <- matrix(as.numeric(object@results$mPsi),
                  nrow=length(object@results$mPsi),
                  ncol=ncol(B), byrow=FALSE)
    M <- solve((diag(ncol(B)) + (t(B) %*% (B/mpi))), 
               (t(B/mpi)))
    af <- t(M %*% useMat)
    
    return(af)
  }
)

