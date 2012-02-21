## METHOD TO READ IN SUMMARY OUTPUT FROM bfrm FOR DIFFERENT MODEL TYPES
#####


setMethod(
  f = ".readResult",
  signature = c("bfrmModel", "character"),
  definition = function(object, outLoc){
    
    theseFiles <- list(sub(".txt", "", list.files(outLoc), fixed=T))
    
    outList <- lapply(theseFiles, fuction(x){
      read.delim(file.path(outLoc, paste(x, ".txt", sep="")), header=F, as.is=T)
    })
    
    newObj <- new("bfrmResult",
                  bfrmModel = object,
                  output = outList)
    
    return(newObj)
  }
)


#####
## SET A SHOW METHOD FOR GENERIC bfrmModel
#####
setMethod(
  f = "show",
  signature = "bfrmResult",
  definition = function(object){
    cat('An object of class "', class(object), '"\n\n', sep="")
    
    these <- slotNames(object)
    cat("Contains slots (class)\n")
    cat("----------------------\n")
    for(this in these)
      cat("  ", this, " (", class(slot(object, this)), ")\n", sep="")
  }
)

