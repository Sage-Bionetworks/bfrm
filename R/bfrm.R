setMethod(
  f = "bfrm",
  signature = "formula",
  definition = function(formula, ...){
    
    Call <- match.call()
    
    y <- eval(Call$formula[[2]])
    x <- as.matrix(eval(Call$formula[[3]]))
    
    args <- list(...)
    
    if(any(names(args) == ""))
      stop("Optional arguments passed for bfrmSetup must be named")
    
    paramSpec <- new("bfrmParam")
    if( length(args) != 0L ){
      for( i in names(args) ){
        slot(paramSpec, i) <- args[[i]]
      }
    }
    
    
    
    if( !(class(y) %in% c("factor", "numeric", "Surv")) )
      stop("Response must either be a numeric vector, a factor, or a Surv object")

    if( class(y) == "Surv" ){
      myObj <- new("bfrmSurvivalModel",
                   call = Call,
                   timeToEvent = y[, 1],
                   censor = y[, 2],
                   data = x,
                   paramSpec = paramSpec)
    } else if( class(y) == "factor" ){
      myObj <- new("bfrmCategoricalModel",
                   call = Call,
                   response = y,
                   data = x,
                   paramSpec = paramSpec)
    } else if( class(y) == "numeric"){
      if( all(unique(y) %in% c(0, 1)) ){
        myObj <- new("bfrmBinaryModel",
                     call = Call,
                     response = y,
                     data = x,
                     paramSpec = paramSpec)
      } else{
        myObj <- new("bfrmLinearModel",
                     call = Call,
                     response = y,
                     data = x,
                     paramSpec = paramSpec)
      }
    }
    
    ## RUN THE WORKER THAT WRITES FILES AND RUNS BFRM BINARY
    outSum <- .bfrmWorker(myObj)
    return(outSum)
  }
)

setMethod(
  f = ".bfrmWorker",
  signature = "bfrmModel",
  definition = function(object){
    
    ## PASS ON TO DISPATCH METHOD DIFFERING BY MODEL TYPE TO FILL IN THE REST OF THE PARAMS
    object <- .writeData(object)

    ## UPDATE DEFAULT SETUP WITH USER SPECIFIED VALUES - AND WRITE OUT THE FILE
    paramLoc <- .writeParam(object@paramSpec)
    
    ## SET UP A LOCATION FOR THE OUTPUT FILES TO BE STORED
    curWd <- getwd()
    outLoc <- tempfile(pattern="output", tmpdir=tempdir(), fileext="")
    dir.create(outLoc)
    setwd(outLoc)
    
    #####
    ## RUN bfrm (ALL THAT IS NEEDED IS THE LOCATION OF THE PARAM FILE)
    system(sprintf("%s %s", system.file(sprintf("/exec/%s/bfrm", Sys.getenv("R_ARCH")), package="bfrm"), paramLoc))
    ## NOW THAT bfrm HAS BEEN CALLED, RETURN SUMMARY OF MODEL RUN
    outSum <- .readResult(object, outLoc)
    #####
    
    setwd(curWd)
    
    return(outSum)
  }
)

