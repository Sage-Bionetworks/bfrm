## GENERIC CLASS DEFINITIONS
##
## AUTHOR: BRIAN M. BOT
#####

## CREATE A CLASS THAT CONTAINS ALL THE BFRM PARAMETER INFORMATION
setClass(
  Class = "bfrmParam",
  
  representation = representation(
    NObservations = "numeric",
    NVariables = "numeric",
    NBinaryResponses = "numeric",
    NCategoricalResponses = "numeric",
    NSurvivalResponses = "numeric",
    NContinuousResponses = "numeric",
    NDesignVariables = "numeric",
    NControlVariables = "numeric",
    NLatentFactors = "numeric",
    DataFile = "character",
    HFile = "character",
    ResponseMaskFile = "character",
    XMaskFile = "character",
    ShapeOfB = "numeric",
    NonGaussianFactors = "numeric",
    PriorPsia = "numeric",
    PriorPsib = "numeric",
    PriorSurvivalPsia = "numeric",
    PriorSurvivalPsib = "numeric",
    PriorRhoMean = "numeric",
    PriorRhoN = "numeric",
    PriorPiMean = "numeric",
    PriorPiN = "numeric",
    PriorTauDesigna = "numeric",
    PriorTauDesignb = "numeric",
    PriorTauResponseBinarya = "numeric",
    PriorTauResponseBinaryb = "numeric",
    PriorTauResponseCategoricala = "numeric",
    PriorTauResponseCategoricalb = "numeric",
    PriorTauResponseSurvivala = "numeric",
    PriorTauResponseSurvivalb = "numeric",
    PriorTauResponseContinuousa = "numeric",
    PriorTauResponseContinuousb = "numeric",
    PriorTauLatenta = "numeric",
    PriorTauLatentb = "numeric",
    PriorInterceptMean = "numeric",
    PriorInterceptVar = "numeric",
    PriorContinuousMean = "numeric",
    PriorContinuousVar = "numeric",
    PriorSurvivalMean = "numeric",
    PriorSurvivalVar = "numeric",
    Evol = "numeric",
    EvolVarIn = "numeric",
    EvolIncludeVariableThreshold = "numeric",
    EvolIncludeFactorThreshold = "numeric",
    EvolMiniumVariablesInFactor = "numeric",
    EvolMaximumFactors = "numeric",
    EvolMaximumVariables = "numeric",
    EvolMaximumVariablesPerIteration = "numeric",
    InclusionMethod = "numeric",
    Burnin = "numeric",
    nMCSamples = "numeric",
    PrintIteration = "numeric",
    PriorAlphaa = "numeric",
    PriorAlphab = "numeric",
    EvolVarInFile = "character"),
  
  prototype = prototype(
    NObservations = 0,
    NVariables = 0,
    NBinaryResponses = 0,
    NCategoricalResponses = 0,
    NSurvivalResponses = 0,
    NContinuousResponses = 0,
    NDesignVariables = 1,
    NControlVariables = 0,
    NLatentFactors = 0,
    ShapeOfB = 2,
    NonGaussianFactors = 1,
    PriorPsia = 10,
    PriorPsib = 2,
    PriorSurvivalPsia = 2,
    PriorSurvivalPsib = 0.5,
    PriorRhoMean = 0.001,
    PriorRhoN = 200,
    PriorPiMean = 0.9,
    PriorPiN = 10,
    PriorTauDesigna = 5,
    PriorTauDesignb = 1,
    PriorTauResponseBinarya = 5,
    PriorTauResponseBinaryb = 1,
    PriorTauResponseCategoricala = 5,
    PriorTauResponseCategoricalb = 1,
    PriorTauResponseSurvivala = 5,
    PriorTauResponseSurvivalb = 1,
    PriorTauResponseContinuousa = 5,
    PriorTauResponseContinuousb = 1,
    PriorTauLatenta = 5,
    PriorTauLatentb = 1,
    PriorInterceptMean = 8,
    PriorInterceptVar = 100,
    PriorContinuousMean = 0,
    PriorContinuousVar = 1,
    PriorSurvivalMean = 2,
    PriorSurvivalVar = 1,
    Evol = 0,
    EvolVarIn = 0,
    EvolIncludeVariableThreshold = 0.75,
    EvolIncludeFactorThreshold = 0.75,
    EvolMiniumVariablesInFactor = 5,
    EvolMaximumFactors = 5,
    EvolMaximumVariables = 100,
    EvolMaximumVariablesPerIteration = 5,
    InclusionMethod = 1,
    Burnin = 2000,
    nMCSamples = 20000,
    PrintIteration = 100,
    PriorAlphaa = 1,
    PriorAlphab = 1)
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
    if( !all(sort(unique(object@response)) %in% c(0, 1)) ){
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
    if( !all(sort(unique(object@censor)) %in% c(0, 1)) ){
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
    response = "numeric")
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
## MODEL RESULT CLASS
#####
setClass(
  Class = "bfrmResult",
  
  representation = representation(
    bfrmModel = "bfrmModel",
    output = "list")
)
