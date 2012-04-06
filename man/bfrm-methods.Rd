\name{bfrm-methods}
\docType{methods}
\alias{bfrm}
\alias{bfrm-methods}
\alias{bfrm,matrix-method}

\title{ Bayesian Factor Regression Modeling (bfrm) }

\description{
  Runs an executable program released by West et al and available here ().  Writes out temporary files, passes them to an executable, and reads summary files back into R.  Temporary files are discarded once the R session is terminated unless specified otherwise.
}

\section{Methods}{
  \describe{
    \item{\code{signature(data = "matrix", ...)}}{}
    \item{\code{data}}{
      numeric matrix with one row per predictor and one column per observation / sample / patient to be passed to the bfrm algorithm
    }
    \item{\code{...}}{
      optional arguments, including:
    }
    \item{\code{design}}{
      design argument(s) to be passed to the bfrm algorithm
    }
    \item{\code{control}}{
      control argument(s) to be passed to the bfrm algorithm (i.e. "assay artifact" variables)
    }
    \item{\code{ybinary}}{
      any binary responses to be considered in the bfrm algorithm
    }
    \item{\code{ycategorical}}{
      any categorical responses to be considered in the bfrm algorithm
    }
    \item{\code{ycontinuous}}{
      any continuous responses to be considered in the bfrm algorithm
    }
    \item{\code{ysurvival}}{
      any time to event responses to be considered in the bfrm algorithm
    }
    \item{\code{censor}}{
      censoring vector to be used with a passed time to event (ysurvival) argument - 0 for censor - 1 for event
    }
    \item{\code{outputDir}}{
      directory for output text files to be stored - default is withing a temporary directory which will be deleted once the R session is terminated.
    }
    \item{\code{other}}{
      named arguments overwriting of defaults specified in the slots of class \code{bfrmParam} - only for advanced users.
    }
  }
}

\section{Value}{
  \describe{
    Return value is of class \code{bfrmResult}
  }
}

\seealso{
  \describe{
    \item{model classes}{ \code{bfrmModel}, \code{evolveModel} }
    \item{methods}{ \code{evolve}, \code{projection} }
    \item{model results}{ \code{bfrmResult} }
  }
}

\keyword{methods}
