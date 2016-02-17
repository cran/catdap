\name{Barplot2WayTable}
\alias{Barplot2WayTable}

\title{Bar Plots for Two-Way Tables}
\description{
  Create bar plots for output two-way tables of catdap1() or catdap2().
}
\usage{
Barplot2WayTable(vname, resvar, exvar=NULL, tway.table, interval=NULL)
}
\arguments{
  \item{vname}{variable names.}
  \item{resvar}{names of the response variables.}
  \item{exvar}{names of the explanatory variables. Default is all variables except \code{resvar}.}
  \item{tway.table}{output \code{tway.table} of 'catdap1()' or 'catdap2()'.}
  \item{interval}{output \code{interval} of 'catdap2()'.}
}
\examples{
## catdap1c()
# Titanic data
resvar <- "Survived"
z1 <- catdap1c(Titanic, resvar)

vname <- names(dimnames(Titanic))
Barplot2WayTable(vname, resvar, , z1$tway.table)

## catdap2()
# Edgar Anderson's Iris Data
# "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width" "Species" 
data(iris)
resvar <- "Petal.Width"
vname <- names(iris)
z2 <- catdap2(iris, c(0, 0, 0, -7, 2), resvar, c(0.1, 0.1, 0.1, 0.1, 0))

exvar <- c("Sepal.Length", "Petal.Length")
Barplot2WayTable(vname, resvar, exvar, z2$tway.table, z2$interval)
}

\keyword{ts}
