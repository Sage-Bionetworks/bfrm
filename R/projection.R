## PROJECT FACTORS ONTO ANOTHER DATASET
## SUMMARIZATION OF RESULTS AS PER MATLAB CODE SENT FROM JOE LUCAS
#####

setMethod(
  f = "projection",
  signature = c("bfrmResult", "matrix"),
  definition = function(factors, newdata){
    
    ## CHECK TO SEE IF NECESSARY OUTPUT FILES ARE PRESENT IN RESULT factors
    if( !all(c("mVariablesIn", "mPostPib", "mA", "mPsi") %in% names(factors@results)) ){
      stop("Required")
    }
    
    ## CHECK IF EITHER SETS OF ROWNAMES ARE NULL - IF SO, MUST ASSUME ORDERED SAME
    odr <- is.null(rownames(factors@model@data))
    ndr <- is.null(rownames(newdata))
    if(odr | ndr){
      warning("rows not named in newdata or original data - must be in same order for projections to be accurate")
      pidNew <- 1:nrow(newdata)
      pidOrig <- 1:nrow(factors@model@data)
    } else{
      pidNew <- rownames(newdata)
      pidOrig <- rownames(factors@model@data)
    }
    
    ## SUBSET VARIABLES IF NECESSARY
    if( any(names(factors@results) == "mVariablesIn") ){
      pidOrig <- pidOrig[as.numeric(factors@results$mVariablesIn)]
      pidNew <- pidNew[as.numeric(factors@results$mVariablesIn)]
    }
    
    if( ! all(pidOrig %in% pidNew) ){
      stop("newdata must contain all features present in original data")
    }
    
    useMat <- newdata[pidOrig, ]
    
    B <- factors@results$mA * factors@results$mPostPib
    mpi <- matrix(as.numeric(factors@results$mPsi),
                  nrow=length(factors@results$mPsi),
                  ncol=ncol(B), byrow=FALSE)
    M <- solve((diag(ncol(B)) + (t(B) %*% (B/mpi))), 
               (t(B/mpi)))
    af <- t(M %*% useMat)
    
    return(af)
  }
)

