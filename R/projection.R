## PROJECT FACTORS ONTO ANOTHER DATASET
## SUMMARIZATION OF RESULTS AS PER MATLAB CODE SENT FROM JOE LUCAS
#####

setMethod(
  f = "projection",
  signature = c("bfrmResult", "matrix"),
  definition = function(object, newdata){
    
    ## CHECK IF EITHER SETS OF ROWNAMES ARE NULL - IF SO, MUST ASSUME ORDERED SAME
    odr <- is.null(rownames(object@bfrmModel@data))
    ndr <- is.null(rownames(newdata))
    if(odr | ndr){
      warning("rows not named in newdata or original data - must be in same order for projections to be accurate")
      pidNew <- 1:nrow(newdata)
      pidOrig <- 1:nrow(object@bfrmModel@data)
    } else{
      pidNew <- rownames(newdata)
      pidOrig <- rownames(object@bfrmModel@data)
    }
    
    ## SUBSET VARIABLES IF NECESSARY
    if( any(names(object@bfrmOutput) == "mVariablesIn") ){
      pidOrig <- pidOrig[as.numeric(object@bfrmOutput$mVariablesIn)]
      pidNew <- pidNew[as.numeric(object@bfrmOutput$mVariablesIn)]
    }
    
    if( ! all(pidOrig %in% pidNew) ){
      stop("newdata must contain all features present in original data")
    }
    
    useMat <- newdata[pidOrig, ]
    
    B <- object@bfrmOutput$mA * object@bfrmOutput$mPostPib
    mpi <- matrix(as.numeric(object@bfrmOutput$mPsi),
                  nrow=length(object@bfrmOutput$mPsi),
                  ncol=ncol(B), byrow=FALSE)
    M <- solve((diag(ncol(B)) + (t(B) %*% (B/mpi))), 
               (t(B/mpi)))
    af <- t(M %*% useMat)
    
    return(af)
  }
)

