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
    .bfrmPlatform(paramLoc)
    ## NOW THAT bfrm HAS BEEN CALLED, RETURN SUMMARY OF MODEL RUN
#    outSum <- .readResult(object, outLoc)
    outSum <- NULL
    #####
    
    setwd(curWd)
    
    return(outSum)
  }
)


#' bfrm platform specific triage
#'
#' Determines if there is a platform specific executable to be run for the user
#'
#' @param paramLoc the location on the local file system of the setup file necessary for
#'   passing to bfrm executable.
setMethod(
  f = ".bfrmPlatform",
  signature = "character",
  definition = function(paramLoc){
    
    ## DEFINE WHICH BINARY TO RUN
    ## myPlat IS ONE OF source, mac.binary.leopard, win.binary
    myPlat <- .Platform$pkgType
    ## FOR DIFFERENTIATION BETWEEN 32 AND 64 BIT OS
    myArch <- .Platform$r_arch
    mySwitch <- paste(myPlat, myArch, sep="")
    
    switch(mySwitch,
           mac.binary.leopardx86_64 = .mac64bfrm(paramLoc),
#           mac.binary.leopardi386   = .mac32bfrm(paramLoc), ## DOES THIS EXIST?
           win.binaryx86_64         = .win64bfrm(paramLoc),
           win.binaryi386           = .win32bfrm(paramLoc),
           source                   = .source64bfrm(paramLoc),
#           sourcex86_64             = .source64bfrm(paramLoc),
#           sourcei386               = .source32bfrm(paramLoc),
           stop("No bfrm executable available for this platform\nSee http://www.stat.duke.edu/research/software/west/bfrm/ for info on available platforms."))
  }
)


#' Platform specific bfrm call
#'
#' Mac-specific executable
#'
#' @param paramLoc the location on the local file system of the setup file necessary for
#'   passing to bfrm executable.
setMethod(
  f = ".mac64bfrm",
  signature = c("character"),
  definition = function(paramLoc){
    pathToExec <- file.path(path.package("bfrm"), "exec")
    system(paste(file.path(pathToExec, "bfrm"), paramLoc, sep=" "))
  }
)

#' Platform specific bfrm call
#'
#' Windows 64-bit specific executable
#'
#' @param paramLoc the location on the local file system of the setup file necessary for
#'   passing to bfrm executable.
setMethod(
  f = ".win64bfrm",
  signature = c("character"),
  definition = function(paramLoc){
    pathToExec <- file.path(path.package("bfrm"), "exec")
    system(paste(file.path(pathToExec, "bfrm64.exe"), paramLoc, sep=" "))
  }
)

#' Platform specific bfrm call
#'
#' Windows 32-bit specific executable
#'
#' @param paramLoc the location on the local file system of the setup file necessary for
#'   passing to bfrm executable.
setMethod(
  f = ".win32bfrm",
  signature = c("character"),
  definition = function(paramLoc){
    pathToExec <- file.path(path.package("bfrm"), "exec")
    system(paste(file.path(pathToExec, "bfrm.exe"), paramLoc, sep=" "))
  }
)

#' Platform specific bfrm call
#'
#' Unix 64-bit specific executable
#'
#' @param paramLoc the location on the local file system of the setup file necessary for
#'   passing to bfrm executable.
setMethod(
  f = ".source64bfrm",
  signature = c("character"),
  definition = function(paramLoc){
    pathToExec <- file.path(path.package("bfrm"), "exec")
    system(paste(file.path(pathToExec, "bfrm64"), paramLoc, sep=" "))
  }
)

#' Platform specific bfrm call
#'
#' Unix 32-bit specific executable
#'
#' @param paramLoc the location on the local file system of the setup file necessary for
#'   passing to bfrm executable.
setMethod(
  f = ".source32bfrm",
  signature = c("character"),
  definition = function(paramLoc){
    pathToExec <- file.path(path.package("bfrm"), "exec")
    system(paste(file.path(pathToExec, "bfrm32"), paramLoc, sep=" "))
  }
)


