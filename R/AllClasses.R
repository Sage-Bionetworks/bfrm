## GENERIC CLASS DEFINITIONS
##
## AUTHOR: BRIAN M. BOT
#####

## CREATE A CLASS THAT CONTAINS ALL THE BFRM PARAMETER INFORMATION
setClass(
  Class = "bfrmParam",
  
  representation = representation(
    nobservations = "numeric",
    nvariables = "numeric",
    nbinaryresponses = "numeric",
    ncategoricalresponses = "numeric",
    nsurvivalresponses = "numeric",
    ncontinuousresponses = "numeric",
    ndesignvariables = "numeric",
    ncontrolvariables = "numeric",
    nlatentfactors = "numeric",
    datafile = "character",
    hfile = "character",
    responsemaskfile = "character",
    xmaskfile = "character",
    shapeofb = "numeric",
    nongaussianfactors = "numeric",
    priorpsia = "numeric",
    priorpsib = "numeric",
    priorsurvivalpsia = "numeric",
    priorsurvivalpsib = "numeric",
    priorrhomean = "numeric",
    priorrhon = "numeric",
    priorpimean = "numeric",
    priorpin = "numeric",
    priortaudesigna = "numeric",
    priortaudesignb = "numeric",
    priortauresponsebinarya = "numeric",
    priortauresponsebinaryb = "numeric",
    priortauresponsecategoricala = "numeric",
    priortauresponsecategoricalb = "numeric",
    priortauresponsesurvivala = "numeric",
    priortauresponsesurvivalb = "numeric",
    priortauresponsecontinuousa = "numeric",
    priortauresponsecontinuousb = "numeric",
    priortaulatenta = "numeric",
    priortaulatentb = "numeric",
    priorinterceptmean = "numeric",
    priorinterceptvar = "numeric",
    priorcontinuousmean = "numeric",
    priorcontinuousvar = "numeric",
    priorsurvivalmean = "numeric",
    priorsurvivalvar = "numeric",
    evol = "numeric",
    evolvarin = "numeric",
    evolincludevariablethreshold = "numeric",
    evolincludefactorthreshold = "numeric",
    evolminiumvariablesinfactor = "numeric",
    evolmaximumfactors = "numeric",
    evolmaximumvariables = "numeric",
    evolmaximumvariablesperiteration = "numeric",
    inclusionmethod = "numeric",
    burnin = "numeric",
    nmcsamples = "numeric",
    printiteration = "numeric",
    prioralphaa = "numeric",
    prioralphab = "numeric",
    evolvarinfile = "character"),
  
  prototype = prototype(
    nobservations = 0,
    nvariables = 0,
    nbinaryresponses = 0,
    ncategoricalresponses = 0,
    nsurvivalresponses = 0,
    ncontinuousresponses = 0,
    ndesignvariables = 1,
    ncontrolvariables = 0,
    nlatentfactors = 0,
    shapeofb = 2,
    nongaussianfactors = 1,
    priorpsia = 10,
    priorpsib = 2,
    priorsurvivalpsia = 2,
    priorsurvivalpsib = 0.5,
    priorrhomean = 0.001,
    priorrhon = 200,
    priorpimean = 0.9,
    priorpin = 10,
    priortaudesigna = 5,
    priortaudesignb = 1,
    priortauresponsebinarya = 5,
    priortauresponsebinaryb = 1,
    priortauresponsecategoricala = 5,
    priortauresponsecategoricalb = 1,
    priortauresponsesurvivala = 5,
    priortauresponsesurvivalb = 1,
    priortauresponsecontinuousa = 5,
    priortauresponsecontinuousb = 1,
    priortaulatenta = 5,
    priortaulatentb = 1,
    priorinterceptmean = 8,
    priorinterceptvar = 100,
    priorcontinuousmean = 0,
    priorcontinuousvar = 1,
    priorsurvivalmean = 2,
    priorsurvivalvar = 1,
    evol = 0,
    evolvarin = 0,
    evolincludevariablethreshold = 0.75,
    evolincludefactorthreshold = 0.75,
    evolminiumvariablesinfactor = 5,
    evolmaximumfactors = 5,
    evolmaximumvariables = 100,
    evolmaximumvariablesperiteration = 5,
    inclusionmethod = 1,
    burnin = 2000,
    nmcsamples = 20000,
    printiteration = 100,
    prioralphaa = 1,
    prioralphab = 1)
  )

## VIRTUAL CLASS THAT WILL BE EXTENDED BY EACH MODEL TYPE


#' ssModel virtual class and the classes that contain it
#'
#' @alias bfrmLinearModel,bfrmBinaryModel,bfrmSurvivalModel
#' @exportClass
setClass(
  Class = "bfrmModel",
  
  representation = representation(
    "VIRTUAL",
    call = "call",
    data = "matrix",
    paramSpec = "bfrmParam")  
)

#####
## INDIVIDUAL MODEL TYPES EXTENDING bfrmModel
#####

## LINEAR MODEL
setClass(
  Class = "bfrmLinearModel",
  contains = "bfrmModel",
  
  representation = representation(
    response = "numeric")
)
setValidity(
  "bfrmLinearModel",
  function(object){
    
    if( ncol(object@data) != length(object@response) ){
      return("number of columns in data does not match number of values in response")
    }
    if( any(is.na(object@response)) ){
      return("NAs not allowed in response vector")
    }
    
    ## IF PASS ABOVE CHECKS THEN RETURN TRUE
    return(TRUE)
  }
)

## BINARY MODEL
setClass(
  Class = "bfrmBinaryModel",
  contains = "bfrmModel",
  
  representation = representation(
    response = "numeric")
)
setValidity(
  "bfrmBinaryModel",
  function(object){
    
    if( ncol(object@data) != length(object@response) ){
      return("number of columns in data does not match number of values in response")
    }
    if( any(is.na(object@response)) ){
      return("NAs not allowed in response vector")
    }
    if( any(!(object@response %in% c(0, 1))) ){
      stop("response for binary must only contain values of 0 and 1")
    }
    
    ## IF PASS ABOVE CHECKS THEN RETURN TRUE
    return(TRUE)
  }
)

## SURVIVAL MODEL
setClass(
  Class = "bfrmSurvivalModel",
  contains = "bfrmModel",
  
  representation = representation(
    timeToEvent = "numeric",
    censor = "numeric")  
)
setValidity(
  "bfrmSurvivalModel",
  function(object){
    
    if( ncol(object@data) != length(object@censor) ){
      return("number of columns in data does not match number of values in censor")
    }
    if( ncol(object@data) != length(object@timeToEvent) ){
      return("number of columns in data does not match number of values in timeToEvent")
    }
    if( any(is.na(object@censor)) ){
      return("NAs not allowed in censor vector")
    }
    if( any(!(object@censor %in% c(0, 1))) ){
      return("censor must only contain values of 0 (censor) and 1 (event)")
    }
    
    ## IF PASS ABOVE CHECKS THEN RETURN TRUE
    return(TRUE)
  }
)

## CATEGORICAL MODEL - NOT AS SURE ABOUT THIS ONE
setClass(
  Class = "bfrmCategoricalModel",
  contains = "bfrmModel",
  
  representation = representation(
    response = "factor")
)
setValidity(
  "bfrmCategoricalModel",
  function(object){
    
    if( ncol(object@data) != length(object@response) ){
      return("number of columns in data does not match number of values in response")
    }
    if( any(is.na(object@response)) ){
      return("NAs not allowed in response vector")
    }
    
    ## IF PASS ABOVE CHECKS THEN RETURN TRUE
    return(TRUE)
  }
)


#####
## MODEL RESULT CLASSES
#####

## VIRTUAL CLASS THAT WILL BE EXTENDED BY EACH MODEL TYPE
setClass(
  Class = "bfrmResult",
  
  representation = representation(
    "VIRTUAL",
    bfrmModel = "bfrmModel",
    bfrmOutput = "list")
)

setClass(
  Class = "bfrmLinearResult",
  contains = "bfrmResult",
  
  representation = representation(
    bfrmModel = "bfrmLinearModel")
)

setClass(
  Class = "bfrmBinaryResult",
  contains = "bfrmResult",
  
  representation = representation(
    bfrmModel = "bfrmBinaryModel")
)

setClass(
  Class = "bfrmSurvivalResult",
  contains = "bfrmResult",
  
  representation = representation(
    bfrmModel = "bfrmSurvivalModel")
)
