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
    hfile = "NA",
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
    nmcsamples = 5000,
    printiteration = 100,
    prioralphaa = 1,
    prioralphab = 1,
    evolvarinfile = "NA")
)

## MODEL CLASS TO DISPATCH ON LATER
setClass(
  Class = "bfrmModel",
  
  representation = representation(
    data = "matrix",
    y = "matrix",
    design = "matrix",
    control = "matrix",
    ymask = "matrix",
    paramSpec = "bfrmParam")  
)
setValidity(
  "bfrmModel",
  function(object){
    
    if( ncol(object@data) != ncol(object@y) ){
      return("number of columns in data does not match number of values in y")
    }
    if( ncol(object@data) != ncol(object@design) ){
      return("number of columns in data does not match number of values in design")
    }
    if( ncol(object@data) != ncol(object@control) ){
      return("number of columns in data does not match number of values in control")
    }
    if( ncol(object@data) != ncol(object@ymask) ){
      return("number of columns in data does not match number of values in censor")
    }
    
    ## IF PASS ABOVE CHECKS THEN RETURN TRUE
    return(TRUE)
  }
)

setClass(
  Class = "evolveModel",
  contains = "bfrmModel",
  
  representation = representation(
    evolveSpecs = "list")
)

#####
## MODEL RESULT CLASS
#####
setClass(
  Class = "bfrmResult",
  
  representation = representation(
    model = "bfrmModel",
    results = "list")
)


#####
## HERE ARE ALL OF THE SHOW METHODS
#####

## SET A SHOW METHOD FOR GENERIC bfrmModel
setMethod(
  f = "show",
  signature = "bfrmResult",
  definition = function(object){
    
#     cat("Call: ", deparse(object@bfrmModel@call), "\n", sep="")
#     cat("Number of features searched : ", ncol(object@bfrmModel@data), "\n", sep="")
#     cat("Number of training samples  : ", sum(object@bfrmModel@training==1), "\n", sep="")
#     if( any(object@bfrmModel@training==0) ){
#       cat("Number of testing samples   : ", sum(object@bfrmModel@training==0), "\n\n", sep="")
#       cat("To access the predictions on the held-out testing dataset, call:\n")
#       cat("  predict(object)\n", sep="")
#       
#       ## ADD OTHER INFO ABOUT THE PREDICTIONS
#     } else{
#       cat("To test this predictive model against a validation set, pass a new feature matrix to:\n")
#       cat("  predict(object, newdata=newFeatureMatrix)\n", sep="")
#     }
    
    these <- slotNames(object)
    cat("\n----------------------\n")
    cat("Contains slots (class)\n")
    cat("----------------------\n")
    for(this in these){
      cat("  ", this, " (", class(slot(object, this)), ")\n", sep="")
      if( class(slot(object,this)) == "list" ){
        theseL <- names(slot(object, this))
        for(thisL in theseL)
          cat("      ", thisL, "\n", sep="")
      }
    }
    
  }
)


## SET A SHOW METHOD FOR GENERIC bfrmModel
setMethod(
  f = "show",
  signature = "bfrmModel",
  definition = function(object){
    cat('An object of class "', class(object), '"\n\n', sep="")
    
    these <- slotNames(object)
    cat("----------------------\n")
    cat("Contains slots (class)\n")
    cat("----------------------\n")
    for(this in these){
      cat("  ", this, " (", class(slot(object, this)), ")\n", sep="")
      if( class(this) == "list" ){
        theseL <- names(object)
        for(thisL in theseL)
          cat("    ", thisL, " (", class(slot(object, this)[[thisL]]), ")\n", sep="")
      }
    }
  }
)

## SET A SHOW METHOD FOR GENERIC bfrmParam
setMethod(
  f = "show",
  signature = "bfrmParam",
  definition = function(object){
    cat('An object of class "', class(object), '"\n\n', sep="")
    
    these <- slotNames(object)
    cat("---------------------------------\n")
    cat("Contains parameter slots (values)\n")
    cat("---------------------------------\n")
    for(this in these){
      cat("  ", this, " (", slot(object, this), ")\n", sep="")
    }
  }
)


