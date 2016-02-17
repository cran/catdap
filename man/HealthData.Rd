\name{HealthData}
\alias{HealthData}

\title{Health Data}
\description{
  Medical data containing both continuous and categorical explanatory variables.
}
\usage{
data(HealthData)
}
\format{
  A data frame with 52 observations on the following 8 variables.

  A part of the source data was recoded according to an input example of original program CATDAP-02.
  In addition we converted 1 into "A" and 2 into "B" of symptoms data, and converted cholesterol data less than 198 into "low" and the others into "high".

  \tabular{rll}{
  [, 1] \tab opthalmo.\tab 1, 2\cr
  [, 2] \tab ecg\tab 1, 2\cr
  [, 3] \tab symptoms\tab A, B\cr
  [, 4] \tab age\tab 49-59\cr
  [, 5] \tab max.press\tab 98-216\cr
  [, 6] \tab min.press\tab 56-120\cr
  [, 7] \tab aortic.wav\tab 6.3-10.2\cr
  [, 8] \tab cholesterol\tab low, high
  }
}
\source{
  Y.Sakamoto, M.Ishiguro and G.Kitagawa (1982)
  \emph{Computer Science Monograph, No.14, CATDAP, A CATEGORICAL DATA ANALYSIS PROGRAM PACKAGE, DATA No.2.}
  The Institute of Statistical Mathematics.

  Y.Sakamoto (1985) \emph{Categorical Data Analysis by AIC, p. 74.} Kluwer Academic publishers.
}

\keyword{datasets}
