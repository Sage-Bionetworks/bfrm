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
    \item{\code{nobservations}:}{Object of class \code{"numeric"} ~~ total number of samples (observations) in data -- SET AUTOMATICALLY BY \code{bfrm} or \code{evolve} FUNCTION CALL}
    \item{\code{nvariables}:}{Object of class \code{"numeric"} ~~ total number of varialbes (features) in data -- SET AUTOMATICALLY BY \code{bfrm} or \code{evolve} FUNCTION CALL}
    \item{\code{nbinaryresponses}:}{Object of class \code{"numeric"} ~~ total number of binary response variables -- SET AUTOMATICALLY BY \code{bfrm} or \code{evolve} FUNCTION CALL}
    \item{\code{ncategoricalresponses}:}{Object of class \code{"numeric"} ~~ total number of categorical response variables -- SET AUTOMATICALLY BY \code{bfrm} or \code{evolve} FUNCTION CALL}
    \item{\code{nsurvivalresponses}:}{Object of class \code{"numeric"} ~~ total number of survival response variables -- SET AUTOMATICALLY BY \code{bfrm} or \code{evolve} FUNCTION CALL}
    \item{\code{ncontinuousresponses}:}{Object of class \code{"numeric"} ~~ total number of continuous response variables -- SET AUTOMATICALLY BY \code{bfrm} or \code{evolve} FUNCTION CALL}
    \item{\code{ndesignvariables}:}{Object of class \code{"numeric"} ~~ number of design covariates, including the intercept -- SET AUTOMATICALLY BY \code{bfrm} or \code{evolve} FUNCTION CALL -- default value is 1 if no design variables are passed}
    \item{\code{ncontrolvariables}:}{Object of class \code{"numeric"} ~~ total number of control, or "assay-artifact", covariates -- SET AUTOMATICALLY BY \code{bfrm} or \code{evolve} FUNCTION CALL -- default value is 0 if no control varialbes are passed}
    \item{\code{nlatentfactors}:}{Object of class \code{"numeric"} ~~ either (1) the number of latent factors in the model for static \code{bfrm} model; or (2) the starting number of latent factors in the model for evolutionary bfrm (\code{evolve}) model -- SET AUTOMATICALLY BY \code{bfrm} or \code{evolve} FUNCTION CALL}
    \item{\code{datafile}:}{Object of class \code{"character"} ~~ name of the datafile passed to the bfrm binary -- SET AUTOMATICALLY BY \code{bfrm} or \code{evolve} FUNCTION CALL}
    \item{\code{hfile}:}{Object of class \code{"character"} ~~ name of the file that contains the design and control variables passed to the bfrm binary -- SET AUTOMATICALLY BY \code{bfrm} or \code{evolve} FUNCTION CALL}
    \item{\code{responsemaskfile}:}{Object of class \code{"character"} ~~ name of the file that indicates response masking passed to the bfrm binary -- SET AUTOMATICALLY BY \code{bfrm} or \code{evolve} FUNCTION CALL}
    \item{\code{xmaskfile}:}{Object of class \code{"character"} ~~ name of the file that indicates masked variables in the data matrix passed to the bfrm binary -- SET AUTOMATICALLY BY \code{bfrm} or \code{evolve} FUNCTION CALL}
    \item{\code{shapeofb}:}{Object of class \code{"numeric"} ~~ This parameter defines the constraints placed on the factor loading matrix B. Takes either 0 (no constraint) or 2 (upper triangular of B set to zero) as its value -- default value is 2}
    \item{\code{nongaussianfactors}:}{Object of class \code{"numeric"} ~~ indicator as to whether a gaussian model (0) or Dirichlet Process (1) model is used for latent factors -- default value is 1}
    \item{\code{priorpsia}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the inverse-Gamma(a,b) prior for elements of Psi, the vector of residual varialbes for all X varianbles -- default value is 2 which is set for Affy gene expression data using rma log2}
    \item{\code{priorpsib}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the inverse-Gamma(a,b) prior for elements of Psi, the vector of residual varialbes for all X varianbles -- default value is 0.005 which is set for Affy gene expression data using rma log2}
    \item{\code{priorsurvivalpsia}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the inverse-gamma(a,b) prior for residual variances of an included survival response variable; right censored survival data are modelled as log-normal, linear regressions -- default value is 2}
    \item{\code{priorsurvivalpsib}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the inverse-gamma(a,b) prior for residual variances of an included survival response variable; right censored survival data are modelled as log-normal, linear regressions -- default value is 0.5}
    \item{\code{priorrhomean}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the Beta(PriorRhoMean* PriorRhoN, (1-PriorRhoMean)*PriorRhoN) prior for the sparsity base rate parameters - the elements of the vector Rho -- defaults value is 0.001}
    \item{\code{priorrhon}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the Beta(PriorRhoMean* PriorRhoN, (1-PriorRhoMean)*PriorRhoN) prior for the sparsity base rate parameters - the elements of the vector Rho -- defaults value is 200}
    \item{\code{priorpimean}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the Beta(PriorPiMean* PriorPiN, (1-PriorPiMean)*PriorPiN) prior for the hierachical components of the prior on non-zero inclusion probabilities -- default value is 0.9}
    \item{\code{priorpin}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the Beta(PriorPiMean* PriorPiN, (1-PriorPiMean)*PriorPiN) prior for the hierachical components of the prior on non-zero inclusion probabilities -- default value is 10}
    \item{\code{priortaudesigna}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the inverse-Gamma(a,b) prior for the variances Tau of the design/control factor effects -- default value is 5}
    \item{\code{priortaudesignb}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the inverse-Gamma(a,b) prior for the variances Tau of the design/control factor effects -- default value is 1}
    \item{\code{priortauresponsebinarya}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the inverse-Gamma(a,b) prior for the variances Tau of the binary response factors -- default value is 5}
    \item{\code{priortauresponsebinaryb}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the inverse-Gamma(a,b) prior for the variances Tau of the binary response factors -- default value is 1}
    \item{\code{priortauresponsecategoricala}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the inverse-Gamma(a,b) prior for the variances Tau of the categorical response factors -- default value is 5}
    \item{\code{priortauresponsecategoricalb}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the inverse-Gamma(a,b) prior for the variances Tau of the categorical response factors -- default value is 1}
    \item{\code{priortauresponsesurvivala}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the inverse-Gamma(a,b) prior for the variances Tau in the of the survival response factors -- default value is 5}
    \item{\code{priortauresponsesurvivalb}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the inverse-Gamma(a,b) prior for the variances Tau in the of the survival response factors -- default value is 1}
    \item{\code{priortauresponsecontinuousa}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the inverse-Gamma(a,b) prior for the variances Tau of the continuous response factors -- default value is 5}
    \item{\code{priortauresponsecontinuousb}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the inverse-Gamma(a,b) prior for the variances Tau of the continuous response factors -- default value is 1}
    \item{\code{priortaulatenta}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the inverse-Gamma(a,b) prior of the variances Tau for the latent factors -- default value is 5}
    \item{\code{priortaulatentb}:}{Object of class \code{"numeric"} ~~ hyper-parameter values for the inverse-Gamma(a,b) prior of the variances Tau for the latent factors -- default value is 1}
    \item{\code{priorinterceptmean}:}{Object of class \code{"numeric"} ~~ prior mean for the intercept (baseline level) of X variables -- default value is 8 based on prototype of Affy gene expression data}
    \item{\code{priorinterceptvar}:}{Object of class \code{"numeric"} ~~ prior variance for the intercept (baseline level) of X variables -- default value is 100 based on prototype of Affy gene expression data}
    \item{\code{priorcontinuousmean}:}{Object of class \code{"numeric"} ~~ prior mean for the intercept (baseline) of any continuous response variables -- default value is 0}
    \item{\code{priorcontinuousvar}:}{Object of class \code{"numeric"} ~~ prior variance for the intercept (baseline) of any continuous response variables -- default value is 1}
    \item{\code{priorsurvivalmean}:}{Object of class \code{"numeric"} ~~ prior mean for the intercept (baseline) of any survival response variables -- default value is 2}
    \item{\code{priorsurvivalvar}:}{Object of class \code{"numeric"} ~~ prior mean for the intercept (baseline) of any survival response variables -- default value is 10}
    \item{\code{evol}:}{Object of class \code{"numeric"} ~~ indicator for evolutionary mode -- SET AUTOMATICALLY BY \code{bfrm} (0) or \code{evolve} (1) FUNCTION CALL}
    \item{\code{evolvarin}:}{Object of class \code{"numeric"} ~~ number of variables used to initialize the evolutionary analysis -- SET AUTOMATICALLY BY \code{evolve} FUNCTION CALL}
    \item{\code{evolvarinfile}:}{Object of class \code{"character"} ~~ indices of the variables included in the initializing set -- SET AUTOMATICALLY BY \code{evolve} FUNCTION CALL -- default value is 1 (first variable in X matrix)}
    \item{\code{evolincludevariablethreshold}:}{Object of class \code{"numeric"} ~~ threshold for bringing a new variables into the model -- default value is 0.75}
    \item{\code{evolincludefactorthreshold}:}{Object of class \code{"numeric"} ~~ threshold for adding a new latent factor into the model -- default value is 0.75}
    \item{\code{evolminiumvariablesinfactor}:}{Object of class \code{"numeric"} ~~ minimum number of variables (genes) showing significant association with a factor in order for the factor to be included in the model -- default value is 5}
    \item{\code{evolmaximumfactors}:}{Object of class \code{"numeric"} ~~ maximum number of latent factors that the final model can have -- default value is 5}
    \item{\code{evolmaximumvariables}:}{Object of class \code{"numeric"} ~~ maximum number of variables the final model can have -- default value is 100}
    \item{\code{evolmaximumvariablesperiteration}:}{Object of class \code{"numeric"} ~~ maximum number of variables that can be added to the model at each iteration -- default value is 5}
    \item{\code{evolmaximumvariablesperfactor}:}{Object of class \code{"numeric"} ~~ maximum number of variables that can be weighted on any one factor -- default value is 15}
    \item{\code{inclusionmethod}:}{Object of class \code{"numeric"} ~~ default value is 1}
    \item{\code{burnin}:}{Object of class \code{"numeric"} ~~ number of burn-in iterations in the MCMC -- default value is 2000}
    \item{\code{nmcsamples}:}{Object of class \code{"numeric"} ~~ number of MCMC iterations -- default value is 5000}
    \item{\code{printiteration}:}{Object of class \code{"numeric"} ~~ number defining how often MCMC iterations are printed to screen -- default value is 100}
    \item{\code{prioralphaa}:}{Object of class \code{"numeric"} ~~ prior parameters for the Gamma prior for Alpha -- default value is 1}
    \item{\code{prioralphab}:}{Object of class \code{"numeric"} ~~ prior parameters for the Gamma prior for Alpha -- default value is 1}
  }
}

\seealso{
  \describe{
    \item{model classes}{ \code{bfrmModel}, \code{evolveModel} }
    \item{methods}{ \code{bfrm}, \code{projection} }
    \item{model results}{ \code{bfrmResult} }
  }
}

\examples{
showClass("bfrmParam")
}

\keyword{classes}
