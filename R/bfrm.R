#' Bayesian Factor Regression Modelling (bfrm) main function call
#'
#' This function takes an object of a derived class of bfrmModel depending
#' on the type of model being run (linear, binary, categorical, or survival).  This
#' function populated necessary input information for the compiled binary
#' to return results.
#'
#' @param object previously constructed object of class derived from \code{bfrmModel}
#' @return an object derived from class bfrmResult dependent on the type of
#'   model that was fit (linear, binary, categorical, survival)
#' @export
#' @examples
#' NEED TO FILL IN EXAMPLES ONCE PACKAGE BUILT AND DATA DIRECTORY AVAILABLE
setMethod(
  f = "bfrm",
  signature = "bfrmModel",
  definition = function(object){
    
    ## FILL IN SETUP INFORMATION
    object@paramSpec@NObservations <- ncol(object@data)
    object@paramSpec@NVariables <- nrow(object@data)

    ## PASS ON TO DISPATCH METHOD DIFFERING BY MODEL TYPE TO FILL IN THE REST OF THE PARAMS
    object <- .writeData(object)

    ## UPDATE DEFAULT SETUP WITH USER SPECIFIED VALUES - AND WRITE OUT THE FILE
    paramLoc <- .writeParam(object@paramSpec)
    
    ## SET UP A LOCATION FOR THE OUTPUT FILES TO BE STORED
    curWd <- getwd()
    outLoc <- file.path(tempdir(), "output")
    dir.create(outLoc)
    setwd(outLoc)
    
    #####
    ## RUN bfrm (ALL THAT IS NEEDED IS THE LOCATION OF THE PARAM FILE)
    .bfrmPlatform(paramLoc)
    ## NOW THAT bfrm HAS BEEN CALLED, RETURN SUMMARY OF MODEL RUN
    outSum <- .readResult(object, outLoc)
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
           sourcex86_64             = .source64bfrm(paramLoc),
           sourcei386               = .source32bfrm(paramLoc),
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
    #    pathToExec <- file.path(path.package("bfrm"), "exec")
    pathToExec <- file.path("/Users/brian/workspace/gitRepos/bfrm/inst", "exec")
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
    #    pathToExec <- file.path(path.package("bfrm"), "exec")
    pathToExec <- file.path("/Users/brian/workspace/gitRepos/bfrm/inst", "exec")
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
    #    pathToExec <- file.path(path.package("bfrm"), "exec")
    pathToExec <- file.path("/Users/brian/workspace/gitRepos/bfrm/inst", "exec")
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
    #    pathToExec <- file.path(path.package("bfrm"), "exec")
    pathToExec <- file.path("/Users/brian/workspace/gitRepos/bfrm/inst", "exec")
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
    #    pathToExec <- file.path(path.package("bfrm"), "exec")
    pathToExec <- file.path("/Users/brian/workspace/gitRepos/bfrm/inst", "exec")
    system(paste(file.path(pathToExec, "bfrm32"), paramLoc, sep=" "))
  }
)


