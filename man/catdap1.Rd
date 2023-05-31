\name{catdap1}
\alias{catdap1}
\alias{catdap1c}

\title{Categorical Data Analysis Program Package 01}
\description{
  Calculates the degree of
  association between all the possible pairs of categorical variables. 
}
\usage{
catdap1(cdata, response.names = NULL, plot = 1, gray.shade = FALSE, ask = TRUE)
catdap1c(ctable, response.names = NULL, plot = 1, gray.shade = FALSE, ask = TRUE)
}
\arguments{
  \item{cdata}{categorical data matrix with variable names on the first row.}
  \item{ctable}{cross-tabulation table with a list of variable names.}
  \item{response.names}{variable names of response variables. If \code{NULL}
    (default), all variables are regarded as response variables.}
  \item{plot}{split directions for each level of the mosaic:
    \describe{
      \item{0 :}{no plot,}
      \item{1 :}{horizontal (default),}
      \item{2 :}{alternating directions, beginning with a vertical split.}
    }
  }
  \item{gray.shade}{A logical value indicating whether the gamma-corrected grey
    palette should be used. If \code{FALSE} (default), any color palette is used.} 
  \item{ask}{logical; if \code{TRUE} (default), the user is asked to confirm
    before a new page is started. If \code{FALSE}, each new plot create a new
    page.} 
}
\details{
  This function is an R-function style clone of Sakamoto's CATDAP-01 program
  for categorical data analysis. CATDAP-01 calculates the degree of
  association between all the possible pairs of categorical variables. 

  The degree of association is evaluated by AIC value. See help(catdap2) for
  details about AIC.

  \code{\link{catdap2}} should be used when the best subset and categorization
  of explanatory variables are sought for. Continuous explanatory variables
  could be explanatory variables in case of catdap2.
}
\value{
  \item{tway.table}{two-way tables and ratio.}
  \item{total}{total number of data with corresponding code of variables.}
  \item{aic}{AIC's of explanatory variables for each response variable.}
  \item{aic.order}{list of explanatory variable numbers arranged in ascending
    order of AIC.}
}
\references{
  Y.Sakamoto and H.Akaike (1978) \emph{Analysis of Cross-Classified Data by AIC}.
  Ann. Inst. Statist. Math., 30, pp.185-197.

  K.Katsura and Y.Sakamoto (1980) \emph{Computer Science Monograph, No.14,
  CATDAP, A Categorical Data Analysis Program Package}.
  The Institute of Statistical Mathematics.

  Y.Sakamoto, M.Ishiguro and G.Kitagawa (1983) \emph{Information Statistics}
  Kyoritsu Shuppan Co., Ltd., Tokyo. (in Japanese)

  Y.Sakamoto (1985) \emph{Categorical Data Analysis by AIC.} Kluwer Academic
  publishers.
}
\examples{
## example 1 (The Japanese National Character)
data(JNcharacter)
response <- c("born.again", "difficult", "pleasure", "women.job", "money")
catdap1(JNcharacter, response)

# or, simply  
data(JNcharacter)
catdap1(JNcharacter)

## example 2 (Titanic data)
# A data set with 2201 observations on 4 variables (Class, Sex, Age and Survived)
# cross-tabulating data
catdap1c(Titanic, "Survived")

# individual data
x <- data.frame(Titanic)
y <- data.matrix(x)
n <- dim(y)[1]
nc <- dim(y)[2]
z <- array(, dim = c(nc-1, sum(y[, 5])))
k <- 1
for (i in 1:n)
  if (y[i, nc] != 0) {
    np <- y[i, nc]
    for (j in 1:(nc-1))
      z[j, k:(k+np-1)] <- dimnames(Titanic)[[j]][[y[i, j]]]
    k <- k + np
  }
data <- data.frame(aperm(array(z, dim = c(4,2201)), c(2,1)),
                   stringsAsFactors = TRUE)
names(data) <- names(dimnames(Titanic))
catdap1(data, "Survived")
}

\keyword{category}
