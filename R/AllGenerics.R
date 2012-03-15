## GENERIC METHOD DEFINITIONS
##
## AUTHOR: BRIAN M. BOT
#####

setGeneric(
  name = "bfrm",
  def = function(data, ...){
    standardGeneric("bfrm")
  }
)

setGeneric(
  name = "projection",
  def = function(object, newdata){
    standardGeneric("projection")
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

