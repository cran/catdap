\name{catdap-package}
\alias{catdap-package}
\alias{catdap}
\docType{package}

\title{
  Categorical Data Analysis Program Package
}
\description{
  R functions for categorical data analysis
}
\details{
  This package provides functions for analyzing multivariate data.
  Dependencies of the distribution of specified variable (response variable) to
  other variables (explanatory variables) are derived and evaluated by AIC
  (Akaike Information Criterion).

  Functions \code{\link{catdap1}} and \code{\link{catdap1c}} are for the
  analysis of categorical data. Every variable is specified as the response
  variable in turn and the goodness of other variables as the explanatory
  variables to the specified variable is evaluated by AIC.

  Function \code{\link{catdap2}} can be applied to data where categorical
  variable and numerical variable are mixed. Specifying one variable as the
  response variable, the dependencies of its distribution on sets of other
  variables are investigated. If the response variable is categorical,
  contingency table analysis method is employed. If the response variable is
  numerical, categorizing the response variable by pooling, the problem is
  reduced to the categorical response variable case. This method eventually
  finds the dependency of the histogram of numerical response variable on sets
  of explanatory variables.

  The Fortran source program codes for above functions were published in
  Sakamoto, Ishiguro and Kitagawa (1983), and Frontiers of Times Series Modeling
  3 : Modeling Seasonality & Periodicity ; ISM (2002), respectively.
  
}
\references{
  Y.Sakamoto and H.Akaike (1978) \emph{Analysis of Cross-Classified Data by AIC.}
  Ann. Inst. Statist. Math., 30, pp.185-197.

  K.Katsura and Y.Sakamoto (1980) \emph{Computer Science Monograph, No.14,
  CATDAP, A Categorical Data Analysis Program Package.}
  The Institute of Statistical Mathematics.

  Y.Sakamoto, M.Ishiguro and G.Kitagawa (1983) \emph{Information Statistics}
  Kyoritsu Shuppan Co., Ltd., Tokyo. (in Japanese)

  Y.Sakamoto (1985) \emph{Model Analysis of Categorical Data.} Kyoritsu Shuppan
  Co., Ltd., Tokyo. (in Japanese)

  Y.Sakamoto (1985) \emph{Categorical Data Analysis by AIC.} Kluwer Academic
   publishers.

  \emph{An AIC-based Tool for Data Visualization} (2015),
  \href{http://hdl.handle.net/10787/3614}{NTT DATA Mathematical Systems Inc.}
  (in Japanese)
}

\keyword{package}

