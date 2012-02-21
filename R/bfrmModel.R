## CONSTRUCTOR FOR MODEL CLASSES
#####

#####
## LINEAR MODELS
#####
setMethod(
  f = "LinearModel",
  signature = c("numeric", "data.frame"),
  definition = function(response, data, ...){
    
    args <- list(...)
    if(any(names(args) == ""))
      stop("Optional arguments passed for sssSetup must be named")
    
    new("bfrmLinearModel",
        response = response,
        data = data,
        paramSpec = new("bfrmParam", ...))
  }
)


#####
## BINARY MODEL
#####
setMethod(
  f = "BinaryModel",
  signature = c("numeric", "data.frame"),
  definition = function(response, data, ...){

    args <- list(...)
    if(any(names(args) == ""))
      stop("Optional arguments passed for sssSetup must be named")
    
    new("sssBinaryModel",
        response = response,
        data = data,
        paramSpec = new("bfrmParam", ...))
  }
)


#####
## SURVIVAL MODEL
#####
setMethod(
  f = "SurvivalModel",
  signature = c("numeric", "numeric", "matrix"),
  definition = function(timeToEvent, censor, data, ...){

    args <- list(...)
    if(any(names(args) == ""))
      stop("Optional arguments passed for sssSetup must be named")
    
    new("sssSurvivalModel",
        timeToEvent = timeToEvent,
        censor = censor,
        data = data,
        paramSpec = new("bfrmParam", ...))
  }
)



#####
## SET A SHOW METHOD FOR GENERIC sssModel
#####
setMethod(
  f = "show",
  signature = "bfrmModel",
  definition = function(object){
    cat('An object of class "', class(object), '"\n\n', sep="")
    
    these <- slotNames(object)
    cat("Contains slots (class)\n")
    cat("----------------------\n")
    for(this in these)
      cat("  ", this, " (", class(slot(object, this)), ")\n", sep="")
  }
)

