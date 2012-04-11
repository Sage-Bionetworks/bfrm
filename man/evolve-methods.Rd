\name{evolve-methods}
\docType{methods}
\alias{evolve}
\alias{evolve-methods}
\alias{evolve,matrix-method}

\title{ Evolutionary mode for Bayesian Factor Regression Modeling (bfrm) }

\description{
  Runs an executable program released by West et al and available here ().  Writes out temporary files, passes them to an executable, and reads summary files back into R.  Temporary files are discarded once the R session is terminated unless specified otherwise. This method is specific to what the authors call evolutionary mode of the more general Bayesian Factor Regression Modeling (bfrm).
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
    \item{\code{burnin}}{
      number of burn-in iterations in the MCMC -- default is 2000
    }
    \item{\code{nmcsamples}}{
      number of MCMC iterations -- default is 5000
    }
    \item{\code{init}}{
      the rownames (or indices) of \code{data} that are to be included in the initializing set of evolutionary mode -- default is the first row of \code{data}
    }
    \item{\code{varThreshold}}{
      this parameter sets the threshold for bringing a new variable into the model - in considering whether to add in new variables (genes) at a given evolutionary analysis step, variables are ranked according to their approximate posterior probability of inclusion at that stage - pne of the two elements of the decision to include some of the most highly ranked variables is then a threshold on this posterior inclusion probability – variables with probabilities below that threshold will not be included -- default is 0.75 (acceptable values from 0-1)
    }
    \item{\code{facThreshold}}{
      this parameter sets the threshold for adding a new latent factor into the model - a new latent factor will be added if and only if at least this number of variables (genes) for that factor have posterior probability of association with the factor that exceed this probability threshold -- default is 0.75 (acceptable values from 0-1)
    }
    \item{\code{maxVarIter}}{
      this parameter sets the maximum number of variables (genes) that can be added to the model at each iteration -- default is 5
    }
    \item{\code{minFacVars}}{
      this parameter sets the minimum number of variables (genes) showing significant association with a factor in order for that factor to be included in the model -- default is 5
    }
    \item{\code{maxFacVars}}{
      this parameter sets the maximum number of variables that can be weighted on any one factor in the evolutionary inclusion steps - this allows the user to limit the number of variables brought into the model for each factor and hence to explore more effectively other factor dimensions -- default is 15
    }
    \item{\code{maxFacs}}{
      this parameter sets the maximum number of latent factors that the final model can have -- default is 5
    }
    \item{\code{maxVars}}{
      this parameter sets the maximum number of variables the final model can have -- default is 100
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
    \item{methods}{ \code{bfrm}, \code{projection} }
    \item{model results}{ \code{bfrmResult} }
  }
}

\keyword{methods}

