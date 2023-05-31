\name{Barplot2WayTable}
\alias{Barplot2WayTable}

\title{Bar Plots for Two-Way Tables}
\description{
  Create bar plots for output two-way tables of catdap1() or catdap2().
}
\usage{
Barplot2WayTable(x, exvar = NULL, gray.shade = FALSE)
}
\arguments{
  \item{x}{an output object of \code{"catdap1"} or \code{"catdap2"}.}
  \item{exvar}{names of the explanatory variables. Default is all variables
    except \code{resvar}.}
  \item{gray.shade}{A logical value indicating whether the gamma-corrected grey
    palette should be used. If \code{FALSE} (default), any color palette is used.} 
}
\details{
  For continuous variables, we assume that
  \eqn{b_1, b_2, \dots, b_{m+1}}{b(1), b(2), \dots, b(m+1)} are boundary values
  of \eqn{m} bins. Output value ranges \eqn{r_i}{r(i)} \eqn{(1 \le i \le m)} are
  defined as follows : 

  \deqn{r_i = \left[ \; b_i,\; b_{i+1}\; \right. ) \;\; \mathrm{for} \;1 \le i < m,}{%
        r(i) = [ b(i), b(i+1) ) for 1 <= i < m,}
  \deqn{r_m = \left[ \; b_m,\; b_{m+1}\; \right] .}{%
        r(m) = [ b(m), b(m+1) ] .}
}

\examples{
# catdap1c (Titanic data)
resvar <- "Survived"
z1 <- catdap1c(Titanic, resvar)

Barplot2WayTable(z1)

# catdap2 (Edgar Anderson's Iris Data)
# "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width" "Species" 
data(iris)
resvar <- "Petal.Width"
z2 <- catdap2(iris, c(0, 0, 0, -7, 2), resvar, c(0.1, 0.1, 0.1, 0.1, 0))

exvar <- c("Sepal.Length", "Petal.Length")
Barplot2WayTable(z2, exvar)
}

\keyword{ts}
