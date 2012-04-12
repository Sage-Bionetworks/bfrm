\name{projection-methods}
\docType{methods}
\alias{projection}
\alias{projection-methods}
\alias{projection,bfrmResult,matrix-method}

\title{ Projecting (\code{bfrm}-generated) factors into new data }

\description{
  Method to project a set of factors given my a result of \code{bfrm} into a new dataset
}

\section{Methods}{
  \describe{
    \item{\code{signature(factors = "bfrmResult", newdata = "matrix")}}{}
    \item{\code{factors}}{
      an object of class \code{bfrmResult} produced as output from a call to \code{bfrm}
    }
    \item{\code{newdata}}{
      a numeric matrix that has rownames (feature names) either (1) including all feature names included in data that generated \code{factors}; or (2) ordered in precisely the same order as the data that generated \code{factors}
    }
  }
}

\section{Value}{
  \describe{
    Return value is of a numeric matrix of projections into \code{newdata}
  }
}

\seealso{
  \describe{
    \item{model classes}{ \code{bfrmModel}, \code{evolveModel} }
    \item{methods}{ \code{bfrm}, \code{evolve} }
    \item{model results}{ \code{bfrmResult} }
  }
}

\keyword{methods}
