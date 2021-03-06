## MAIN BFRM CALL AND WORKER
## INTEGRATES WITH COMPILED C++ CODE VIA EXECUTABLE CALL
#####


setMethod(
  f = "bfrm",
  signature = "matrix",
  definition = function(data, ...){
    
    x <- data
    args <- list(...)
    
    #####
    ## ARGUMENT PARSING
    #####
    if(any(names(args) == ""))
      stop("Optional arguments passed for bfrmParam must be named")
    
    if( any(names(args) == "design") ){
      design <- args[["design"]]
      args[["design"]] <- NULL
      if( is.numeric(design) & !is.matrix(design) ){
        design <- matrix(design, ncol=length(design))
      }
    } else{
      design <- matrix(nrow=0, ncol=ncol(x))
    }

    if( any(names(args) == "control") ){
      control <- args[["control"]]
      args[["control"]] <- NULL
      if( is.numeric(control) & !is.matrix(control) ){
        control <- matrix(control, ncol=length(control))
      }
    } else{
      control <- matrix(nrow=0, ncol=ncol(x))
    }
    
    if( any(names(args) == "outputDir") ){
      outputDir <- args[["outputDir"]]
      args[["outputDir"]] <- NULL
    } else{
      outputDir <- tempfile(pattern="output", tmpdir=tempdir(), fileext="")
      dir.create(outputDir)
    }
    
    
    ## SET UP THE PARAMETERS FILE FOR PASS TO C++ EXECUTABLE
    paramSpec <- new("bfrmParam")
    if( length(args) != 0L ){
      for( i in names(args) ){
        slot(paramSpec, i) <- args[[i]]
      }
    }
    
    myObj <- new("bfrmModel",
                 data = x,
                 design = design,
                 control = control,
                 paramSpec = paramSpec)
    
    ## PASS ON TO DISPATCH METHOD DIFFERING BY MODEL TYPE TO FILL IN THE REST OF THE PARAMS
    myObj <- .writeData(myObj)
    
    ## UPDATE DEFAULT SETUP WITH USER SPECIFIED VALUES - AND WRITE OUT THE FILE
    paramLoc <- .writeParam(myObj@paramSpec)
    
    ## SET UP A LOCATION FOR THE OUTPUT FILES TO BE STORED
    curWd <- getwd()
    setwd(outputDir)
    
    #####
    ## RUN bfrm (ALL THAT IS NEEDED IS THE LOCATION OF THE PARAM FILE)
    tmpArch <- Sys.getenv("R_ARCH")
    tmpPlat <- .Platform$pkgType
    if( tmpArch=="" & tmpPlat=="mac.binary.leopard" ){
      tmpArch <- "i386"
    }
    system(sprintf("%s %s", system.file(sprintf("/exec/%s/bfrm", tmpArch), package="bfrm"), paramLoc))
    
    ## NOW THAT bfrm HAS BEEN CALLED, RETURN SUMMARY OF MODEL RUN
    outSum <- .readResult(myObj, outputDir)
    #####
    
    setwd(curWd)
    
    return(outSum)
  }
)

