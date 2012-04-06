\name{bfrmResult-class}
\Rdversion{1.1}
\docType{class}
\alias{bfrmResult-class}
\alias{projection,bfrmResult,matrix-method}

\title{Class \code{"bfrmResult"}}

\description{
  Results as well as record of inputs to the executable bfrm call
}

\section{Objects from the Class}{
  Value returned as a result of a call to either \code{bfrm} or \code{evolve}
}

\section{Slots}{
  \describe{
    \item{\code{model}:}{Object of class \code{"bfrmModel"} - record of inputs to the \code{bfrm} or \code{evolve} call }
    \item{\code{results}:}{Object of class \code{"list"} - list of data matrices output as a result of the executable call }
  }
}

\seealso{
}

\examples{
showClass("bfrmResult")
}

\keyword{classes}
