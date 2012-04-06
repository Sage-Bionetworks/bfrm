\name{bfrmParam-class}
\Rdversion{1.1}
\docType{class}
\alias{bfrmParam-class}

\title{Class \code{"bfrmParam"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("bfrmParam", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{nobservations}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nvariables}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nbinaryresponses}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{ncategoricalresponses}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nsurvivalresponses}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{ncontinuousresponses}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{ndesignvariables}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{ncontrolvariables}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nlatentfactors}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{datafile}:}{Object of class \code{"character"} ~~ }
    \item{\code{hfile}:}{Object of class \code{"character"} ~~ }
    \item{\code{responsemaskfile}:}{Object of class \code{"character"} ~~ }
    \item{\code{xmaskfile}:}{Object of class \code{"character"} ~~ }
    \item{\code{shapeofb}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nongaussianfactors}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priorpsia}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priorpsib}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priorsurvivalpsia}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priorsurvivalpsib}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priorrhomean}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priorrhon}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priorpimean}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priorpin}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priortaudesigna}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priortaudesignb}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priortauresponsebinarya}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priortauresponsebinaryb}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priortauresponsecategoricala}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priortauresponsecategoricalb}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priortauresponsesurvivala}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priortauresponsesurvivalb}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priortauresponsecontinuousa}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priortauresponsecontinuousb}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priortaulatenta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priortaulatentb}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priorinterceptmean}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priorinterceptvar}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priorcontinuousmean}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priorcontinuousvar}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priorsurvivalmean}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priorsurvivalvar}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{evol}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{evolvarin}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{evolincludevariablethreshold}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{evolincludefactorthreshold}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{evolminiumvariablesinfactor}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{evolmaximumfactors}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{evolmaximumvariables}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{evolmaximumvariablesperiteration}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{inclusionmethod}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{burnin}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nmcsamples}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{printiteration}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{prioralphaa}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{prioralphab}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{evolvarinfile}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
No methods defined with class "bfrmParam" in the signature.
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("bfrmParam")
}
\keyword{classes}
