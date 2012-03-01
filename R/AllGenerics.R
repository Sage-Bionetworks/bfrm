## GENERIC METHOD DEFINITIONS
##
## AUTHOR: BRIAN M. BOT
#####

setGeneric(
  name = "bfrm",
  def = function(formula, ...){
    standardGeneric("bfrm")
  }
)
setGeneric(
  name = ".bfrmWorker",
  def = function(object){
    standardGeneric(".bfrmWorker")
  }
)

#####
## bfrm DISPATCH CALLS
#####
setGeneric(
  name = ".bfrmPlatform",
  def = function(paramLoc){
    standardGeneric(".bfrmPlatform")
  }
)


#####
## SPECIFIC bfrm EXECUTABLE CALLS
#####
setGeneric(
  name = ".mac64bfrm",
  def = function(paramLoc){
    standardGeneric(".mac64bfrm")
  }
)
setGeneric(
  name = ".win64bfrm",
  def = function(paramLoc){
    standardGeneric(".win64bfrm")
  }
)
setGeneric(
  name = ".win32bfrm",
  def = function(paramLoc){
    standardGeneric(".win32bfrm")
  }
)
setGeneric(
  name = ".source64bfrm",
  def = function(paramLoc){
    standardGeneric(".source64bfrm")
  }
)
setGeneric(
  name = ".source32bfrm",
  def = function(paramLoc){
    standardGeneric(".source32bfrm")
  }
)



#####
## METHODS FOR WRITING OUT AND READING IN FILES
#####
setGeneric(
  name = ".writeParam",
  def = function(paramSpec){
    standardGeneric(".writeParam")
  }
)
setGeneric(
  name = ".writeData",
  def = function(object){
    standardGeneric(".writeData")
  }
)

setGeneric(
  name = ".readResult",
  def = function(object, outLoc){
    standardGeneric(".readResult")
  }
)




#####
## SET A SHOW METHOD FOR GENERIC bfrmModel
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
