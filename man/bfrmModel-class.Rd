\name{bfrmModel-class}
\Rdversion{1.1}
\docType{class}
\alias{bfrmModel}
\alias{bfrmModel-class}

\title{Class \code{"bfrmModel"}}

\description{
  Class is the created from the call to the \code{bfrm} or \code{evolve} methods and is included as an output in \code{bfrmResult}
}

\section{Objects from the Class}{
  Objects can be created by \code{bfrm} or \code{evolve} methods.
}

\section{Slots}{
  \describe{
    \item{\code{data}:}{Object of class \code{"matrix"} - numeric matrix with one row per predictor and one column per observation / sample / patient }
    \item{\code{design}:}{Object of class \code{"matrix"} - numeric matrix of design covariate(s) as passed to \code{bfrm} method }
    \item{\code{control}:}{Object of class \code{"matrix"} - numeric matrix of control covariate(s) as passed to \code{bfrm} method }
    \item{\code{paramSpec}:}{Object of class \code{"bfrmParam"} an object specifying further optional parameters to be passed to the \code{bfrm} method call (see \code{bfrmParam} for further information) }
  }
}

\seealso{
  \describe{
    \item{model classes}{ \code{evolveModel} }
    \item{methods}{ \code{bfrm}, \code{evolve}, \code{projection} }
    \item{model results}{ \code{bfrmResult} }
  }
}

\examples{
showClass("bfrmModel")
}

\keyword{classes}
