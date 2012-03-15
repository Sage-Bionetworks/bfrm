## MAIN BFRM CALL AND WORKER
## INTEGRATES WITH COMPILED C++ CODE VIA EXECUTABLE CALL
#####


setMethod(
  f = "bfrm",
  signature = "matrix",
  definition = function(data, ...){
    
    x <- data
    args <- list(...)
    
    if(any(names(args) == ""))
      stop("Optional arguments passed for bfrmSetup must be named")
    
    if( any(names(args) == "design") ){
      design <- args[["design"]]
      args[["design"]] <- NULL
    } else{
      design <- matrix(nrow=0, ncol=ncol(x))
    }

    if( any(names(args) == "censor") ){
      censor <- args[["censor"]]
      args[["censor"]] <- NULL
    } else{
      censor <- matrix(nrow=0, ncol=ncol(x))
    }
    
    if( any(names(args) == "outputDir") ){
      outputDir <- args[["outputDir"]]
      args[["outputDir"]] <- NULL
    } else{
      outputDir <- tempfile(pattern="output", tmpdir=tempdir(), fileext="")
      dir.create(outputDir)
    }
    
    ## DETERMINE IF THERE ARE RESPONSE MATRICES PASSES
    ## DETERMINE MISSINGNESS / CENSORING
    #####
    ## HACK IN GOOFY NOTATION FOR HOW BFRM SOURCE CODE EXPECTS CENSORED OBSERVATIONS
    ## 0 IS OBSERVED - 1 IS MISSING (NON-OBSERVED) - AND 2 IS A CENSOR
    y <- matrix(nrow=0, ncol=ncol(x))
    ymask <- matrix(nrow=0, ncol=ncol(x))
    for( i in c("ybinary", "ycategorical", "ysurvival", "ycontinuous") ){
      
      ## LOOK FOR Y VARS
      if( any(names(args) == i) ){
        assign(i, args[[i]])
        args[[i]] <- NULL
        if( is.numeric(get(i)) ){
          assign(i, matrix(get(i), ncol=length(get(i))))
        }
      } else{
        assign(i, matrix(nrow=0, ncol=ncol(x)))
      }
      
      if( !is.matrix(get(i)) ){
        stop(paste(i, "must be numeric vector or matrix"))
      }
      
      if( length(get(i)) > 0 ){
        if( ncol(get(i)) != ncol(data) ){
          stop(paste("number of observations in data matrix not equal to number in", i))
        }
        
        y <- rbind(y, get(i))
        
        if( i == "ysurvival" ){
          ymask <- rbind(ymask, ifelse(censor==0, 2, 0))
        } else{
          ymask <- rbind(ymask, ifelse(is.na(get(i)), 1, 0))
        }
      }
      
    }
    ## END OF SECTION DEFINING RESPONSES
    
    paramSpec <- new("bfrmParam")
    if( length(args) != 0L ){
      for( i in names(args) ){
        slot(paramSpec, i) <- args[[i]]
      }
    }
    
    for( i in c("ybinary", "ycategorical", "ysurvival", "ycontinuous") ){
      slot(paramSpec, paste(sub("y", "n", i), "responses", sep="")) <- nrow(get(i))
    }
    
    myObj <- new("bfrmModel",
                 data = x,
                 y = y,
                 design = design,
                 ymask = ymask,
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
    if(tmpArch=="")
      tmpArch <- "i386"
    system(sprintf("%s %s", system.file(sprintf("/exec/%s/bfrm", tmpArch), package="bfrm"), paramLoc))
    
    ## NOW THAT bfrm HAS BEEN CALLED, RETURN SUMMARY OF MODEL RUN
    outSum <- .readResult(myObj, outputDir)
    #####
    
    setwd(curWd)
    
    return(outSum)
  }
)





#####
## CALLS TO GET BFRM INPUTS INTO CORRECT FORMAT
#####
# setMethod(
#   f = "bfrm",
#   signature = c("data.frame", "ANY", "ANY", "ANY", "ANY"),
#   definition = function(data, ycontinuous, ybinary, ycategorical, ysurvival, ...){
#     bfrm(as.matrix(data), ycontinuous, ybinary, ycategorical, ysurvival, ...)
#   }
# )
# setMethod(
#   f = "bfrm",
#   signature = c("ANY", "numeric", "ANY", "ANY", "ANY"),
#   definition = function(data, ycontinuous, ybinary, ycategorical, ysurvival, ...){
#     bfrm(data, matrix(ycontinuous, ncol=length(ycontinuous)), ybinary, ycategorical, ysurvival, ...)
#   }
# )
# setMethod(
#   f = "bfrm",
#   signature = c("ANY", "ANY", "numeric", "ANY", "ANY"),
#   definition = function(data, ycontinuous, ybinary, ycategorical, ysurvival, ...){
#     bfrm(data, ycontinuous, matrix(ybinary, ncol=length(ybinary)), ycategorical, ysurvival, ...)
#   }
# )
# setMethod(
#   f = "bfrm",
#   signature = c("ANY", "ANY", "ANY", "numeric", "ANY"),
#   definition = function(data, ycontinuous, ybinary, ycategorical, ysurvival, ...){
#     bfrm(data, ycontinuous, ybinary, matrix(ycategorical, ncol=length(ycategorical)), ysurvival, ...)
#   }
# )
# setMethod(
#   f = "bfrm",
#   signature = c("ANY", "ANY", "ANY", "ANY", "numeric"),
#   definition = function(data, ycontinuous, ybinary, ycategorical, ysurvival, ...){
#     bfrm(data, ycontinuous, ybinary, ycategorical, matrix(ysurvival, ncol=length(ysurvival)), ...)
#   }
# )


