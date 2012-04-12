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
  name = "evolve",
  def = function(data, ...){
    standardGeneric("evolve")
  }
)

setGeneric(
  name = "projection",
  def = function(factors, newdata){
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
  def = function(object, outputDir){
    standardGeneric(".readResult")
  }
)

