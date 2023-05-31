\name{MissingHealthData}
\alias{MissingHealthData}

\title{Health Data with Missing Values}
\description{
  Medical data containing both categorical variables and continuous variables,
  the latter include two variables with missing values.
}
\usage{
data(MissingHealthData)
}
\format{
  A data frame with 52 observations on the following 8 variables.

  A part of the source data was recoded according to an input example of
  original program CATDAP-02. In addition, we converted 1 into 'A' and 2 into
  'B' of symptoms data, and converted cholesterol data less than 198 into 'low'
  and the others into 'high'.

  \tabular{rlll}{
  [, 1] \tab opthalmo. \tab 1, 2\cr
  [, 2] \tab ecg \tab\tab 1, 2\cr
  [, 3] \tab symptoms \tab\tab A, B\cr
  [, 4] \tab age \tab\tab 49-59\cr
  [, 5] \tab max.press \tab\tab 98-216, 300 (missing value)\cr
  [, 6] \tab min.press \tab\tab 56-120, 300 (missing value)\cr
  [, 7] \tab aortic.wav \tab\tab 6.3-10.2\cr
  [, 8] \tab cholesterol \tab\tab low, high
  }
}
\source{
  Y.Sakamoto, M.Ishiguro and G.Kitagawa (1980)
  \emph{Computer Science Monograph, No.14, CATDAP, A CATEGORICAL DATA ANALYSIS
  PROGRAM PACKAGE, DATA No.2.} The Institute of Statistical Mathematics.

  Y.Sakamoto (1985) \emph{Categorical Data Analysis by AIC, p. 74.} Kluwer
  Academic publishers.
}

\keyword{datasets}
