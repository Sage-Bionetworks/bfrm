## GENERIC METHOD DEFINITIONS
##
## AUTHOR: BRIAN M. BOT
#####

setGeneric(
  name = "bfrm",
  def = function(object){
    standardGeneric("bfrm")
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
## CLASS CONSTRUCTORS
#####
setGeneric(
  name = "LinearModel",
  def = function(response, data, ...){
    standardGeneric("LinearModel")
  }
)
setGeneric(
  name = "BinaryModel",
  def = function(response, data, ...){
    standardGeneric("BinaryModel")
  }
)
setGeneric(
  name = "CategoricalModel",
  def = function(response, data, ...){
    standardGeneric("CategoricalModel")
  }
)
setGeneric(
  name = "SurvivalModel",
  def = function(timeToEvent, censor, data, ...){
    standardGeneric("SurvivalModel")
  }
)
