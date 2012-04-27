## SPECIAL CASE OF BFRM IN WHICH EVOLUTIONARY MODE IS ENABLED
## INTEGRATES WITH COMPILED C++ CODE VIA EXECUTABLE CALL
#####


setMethod(
  f = "evolve",
  signature = "matrix",
  definition = function(data, ...){
    
    x <- data
    args <- list(...)
    
    #####
    ## PARAMETERS FOR EVOLUTION BY DEFINITION
    #####
    args[["evol"]] <- 1
    args[["nlatentfactors"]] <- 1
    
    #####
    ## ARGUMENT PARSING
    #####
    if(any(names(args) == ""))
      stop("Optional arguments passed for bfrmParam must be named")
    
    if( any(names(args) == "outputDir") ){
      outputDir <- args[["outputDir"]]
      args[["outputDir"]] <- NULL
    } else{
      outputDir <- tempfile(pattern="output", tmpdir=tempdir(), fileext="")
      dir.create(outputDir)
    }
    
    ## EVOLUTIONARY MODE SPECIFIC ARGUMENTS THAT ARE SUGGESTED
    if( any(names(args) == "init") ){
      init <- args[["init"]]
      args[["init"]] <- NULL
      if( length(init) > dim(x)[1] ){
        stop("More init values than rows in data")
      }
      if( class(init) == "character" ){
        if( any(!(init %in% rownames(x))) ){
          stop("init contains character value(s) not included as a rowname in data")
        }
        init <- which(rownames(x) %in% init)
      } else{
        if( !(class(init) %in% c("numeric", "integer")) ){
          stop("init must either be indices or names rows of data")
        } else{
          if( any(init>nrow(x) | init<1) ){
            stop(paste("init values not in the range of 0 to ", nrow(x), sep=""))
          }
        }
      }
    } else{
      init <- 1
    }
    initLoc <- tempfile(pattern="varinfile", tmpdir=tempdir(), fileext=".txt")
    write.table(init, file=initLoc, sep="\t", quote=F, row.names=F, col.names=F)
    args[["evolvarinfile"]] <- initLoc
    args[["evolvarin"]] <- length(init)
    
    if( any(names(args) == "varThreshold") ){
      if( !is.numeric(args[["varThreshold"]]) | length(args[["varThreshold"]]) != 1L ){
        stop("varThreshold must be numeric and between 0 and 1")
      }
      if( args[["varThreshold"]] < 0 | args[["varThreshold"]] > 1 ){
        stop("varThreshold must be numeric and between 0 and 1")
      }
      
      args[["evolincludevariablethreshold"]] <- args[["varThreshold"]]
      args[["varThreshold"]] <- NULL
    }
    if( any(names(args) == "facThreshold") ){
      if( !is.numeric(args[["facThreshold"]]) | length(args[["facThreshold"]]) != 1L ){
        stop("facThreshold must be numeric and between 0 and 1")
      }
      if( args[["facThreshold"]] < 0 | args[["facThreshold"]] > 1 ){
        stop("facThreshold must be numeric and between 0 and 1")
      }
      
      args[["evolincludefactorthreshold"]] <- args[["facThreshold"]]
      args[["facThreshold"]] <- NULL
    }
    
    if( any(names(args) == "maxVarIter") ){
      if( !is.numeric(args[["maxVarIter"]]) | length(args[["maxVarIter"]]) != 1L ){
        stop("maxVarIter must be a single positive integer")
      }
      if( args[["maxVarIter"]] <= 0 ){
        stop("maxVarIter must be a single positive integer")
      }
      
      args[["evolmaximumvariablesperiteration"]] <- args[["maxVarIter"]]
      args[["maxVarIter"]] <- NULL
    }
    
    if( any(names(args) == "maxFacVars") ){
      if( !is.numeric(args[["maxFacVars"]]) | length(args[["maxFacVars"]]) != 1L ){
        stop("maxFacVars must be a single positive integer")
      }
      if( args[["maxFacVars"]] <= 0 ){
        stop("maxFacVars must be a single positive integer")
      }
      
      args[["evolmaximumvariablesperfactor"]] <- args[["maxFacVars"]]
      args[["maxFacVars"]] <- NULL
    }
    
    if( any(names(args) == "minFacVars") ){
      if( !is.numeric(args[["minFacVars"]]) | length(args[["minFacVars"]]) != 1L ){
        stop("minFacVars must be a single non-negative integer")
      }
      if( args[["minFacVars"]] < 0 ){
        stop("minFacVars must be a single non-negative integer")
      }
      
      args[["evolminimumvariablesinfactor"]] <- args[["minFacVars"]]
      args[["minFacVars"]] <- NULL
    }
    
    if( any(names(args) == "maxFacs") ){
      if( !is.numeric(args[["maxFacs"]]) | length(args[["maxFacs"]]) != 1L ){
        stop("maxFacs must be a single positive integer")
      }
      if( args[["maxFacs"]] < 0 ){
        stop("maxFacs must be a single positive integer")
      }
      
      args[["evolmaximumfactors"]] <- args[["maxFacs"]]
      args[["maxFacs"]] <- NULL
    }
    
    if( any(names(args) == "maxVars") ){
      if( !is.numeric(args[["maxVars"]]) | length(args[["maxVars"]]) != 1L ){
        stop("maxVars must be a single positive integer")
      }
      if( args[["maxVars"]] < 0 ){
        stop("maxVars must be a single positive integer")
      }
      
      args[["evolmaximumvariables"]] <- args[["maxVars"]]
      args[["maxVars"]] <- NULL
    }
    
    
    ## MORE GENERIC BFRM ARGUMENTS TO BE PASSED
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
    #####
    ## END OF ARGUMENT PARSING
    #####
    
    ## SET UP THE PARAMETERS FILE FOR PASS TO C++ EXECUTABLE
    paramSpec <- new("bfrmParam")
    slot(paramSpec, "evol") <- 1
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

