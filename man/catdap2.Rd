\name{catdap2}
\alias{catdap2}

\title{Categorical Data Analysis Program Package 02}
\description{
  Search for the best single explanatory variable and detect the best subset of
  explanatory variables.
}
\usage{
catdap2(data, pool = NULL, response.name, accuracy = NULL, nvar = NULL,
        additional.output = NULL, missingmark = NULL, pa1 = 1, pa2 = 4, pa3 = 10,
        print.level = 0, plot = 1)
}
\arguments{
  \item{data}{data matrix with variable names on the first row.}
  \item{pool}{the ways of pooling to categorize each variable must be specified
    by integer parameters:
    \describe{
      \item{(-m) < 0 :}{\eqn{m}-bin histogram is employed to describe the
        distribution of continuous response variable (this option is valid only
        for the response variable),}
      \item{0 :}{equally spaced pooling via a top-down algorithm,}
      \item{1 :}{unequally spaced pooling via a bottom-up algorithm (default),}
      \item{2 :}{no pooling for discrete variables.}
    }
  }
  \item{response.name}{variable name of the response variable.}
  \item{accuracy}{minimum width for the discretization for each variable.}
  \item{nvar}{number of variables to be retained for the analysis of
    multidimensional tables. Default is the number of variables in \code{data}.}
  \item{additional.output}{list of sets of explanatory variable names for
    additional output.}
  \item{missingmark}{positive number for handling missing value. See 'Details'.}
  \item{pa1, pa2, pa3}{control parameter for size of the working area. If error
    message is output, please change the value of parameter according to it.}
  \item{print.level}{this argument determines the level of output printing. The
    default value of '\code{0}' means that lists of "AIC's of the models with k
    explanatory variables (k=1,2,...)" are printed. A value of '\code{1}' means
    that those lists are not printed and "Summary of subsets of explanatory
    variables" within the top 30 is listed.}
  \item{plot}{split directions for each level of the mosaic:
    \describe{
      \item{0 :}{no plot,}
      \item{1 :}{horizontal (default),}
      \item{2 :}{alternating directions, beginning with a vertical split.}
    }
  }
}
\details{
  This function is an R-function style clone of Sakamoto's CATDAP-02 program
  for categorical data analysis. CATDAP-02 can be used to search for the best
  subset of explanatory variables which have the most effective information on a
  specified response variable. Continuous explanatory variables could be
  explanatory variables. In that case CATDAP-02 searches for optimal
  categorization of continuous values.

  The basic statistic adopted is obtained by the application of the
  statistic AIC to the models.

  \eqn{E} denotes the response variable and \eqn{F} denotes candidate
  explanatory variable, and their cell frequencies by
  \eqn{n_E(i) (i\in(E))}{n_E(i) (i in (E))} and
  \eqn{n_F(j) (j\in(F))}{n_F(j) (j in (F))}. The cross frequency is denoted by
  \eqn{n_F,F(i,j)}{n_F,F(i,j)} \eqn{(i,j\in(E,F))}{(i,j in (E,F))}. To measure
  the strength of dependence of a specific set of response variables \eqn{E} on
  the explanatory variable \eqn{F}, we use the following statistic:

  \deqn{AIC(E;F) = -2\sum_{i,j\in(E,F)} n_E,F(i,j)\ ln\{n_E,F(i,j)/(n_F(j)\} + 2(C_E-1)C_F,\ \ (1)}{%
        AIC(E;F) = -2 \sum _{i,j in (E,F)} n_E,F(i,j) ln{n_E,F(i,j)/(n_F(j)} + 2(C_E-1)C_F,  (1)}

  where \eqn{C_E} and \eqn{C_F} denote the total number of categories of the
  corresponding sets of variables, respectively.

  The selection of the best subset of explanatory variables is realized by the
  search for \eqn{F} which gives the minimum \eqn{AIC(E;F)}.

  In case of \eqn{F=\phi}, the formula (1) reduces to

    \deqn{AIC(E;\phi) = -2\sum_{i\in(E)} n_E(i)\ ln\{n_E(i)/n\} + 2(C_E-1).}{%
          AIC(E;\phi) = -2 \sum _{i in (E)} n_E(i) ln{n_E(i)/n} + 2(C_E-1).}

  Here it is assumed that \eqn{C_\phi=1} and \eqn{n_\phi(1)=n}.

  Sakamoto's original CATDAP outputs \eqn{AIC(E;F) - AIC(E;\phi)} as the AIC
  value instead of \eqn{AIC(E;F)}. By this way the positive value of AIC
  indicates that the variable \eqn{F} is judged to be useless as the explanatory
  variable of the \eqn{E}.

  On the other hand, this policy make impossible to compare the goodness of the
  CATDAP model with other models, logit models for example. 

  Considering the convenience of users, present "R version CATDAP" provides not
  only \eqn{AIC = AIC(E;F) - AIC(E;\phi)}, but \eqn{AIC(E;\phi)}, either. The
  latter value is given as base_AIC in the output.

  Users could recover \eqn{AIC(E;F)} by adding AIC and base_AIC.

  \code{missingmark} enables missing value handling.
  When a positive values, say \eqn{1000}, is set here, any value, say \eqn{x},
  greater than or equal to \eqn{1000} is treated as a missing value. If
  \eqn{1000 \le x < 2000}{1000 <= x < 2000}, \eqn{x} is treated as a missing
  value of the 1st type. If \eqn{2000 \le x < 3000}{2000 <= x < 3000}, \eqn{x}
  is treated as a missing value of the 2nd type, and so on. Generally speaking,
  any \eqn{x} that \eqn{1000k \le x < 1000(k+1)}{1000*k <= x < 1000*(k+1)} is
  treated as the \eqn{k}-th type missing value. Users are referred to the
  reference for the technical details of the missing value handling procedure.

  For continuous variables, we assume that
  \eqn{b_1, b_2, \dots, b_{m+1}}{b(1), b(2), \dots, b(m+1)} are boundary values
  of \eqn{m} bins. Output value ranges \eqn{r_i}{r(i)} \eqn{(1 \le i \le m)} are
  defined as follows : 

  \deqn{r_i = \left[ \; b_i,\; b_{i+1}\; \right. ) \;\; \mathrm{for} \;1 \le i < m,}{%
        r(i) = [ b(i), b(i+1) ) for 1 <= i < m,}
  \deqn{r_m = \left[ \; b_m,\; b_{m+1}\; \right] .}{%
        r(m) = [ b(m), b(m+1) ] .}

  Specifically, for continuous response variable \eqn{V},

  \deqn{r_i = \left[ \; x_{min} + (i-1)*s,\; x_{min} + i*s \; \right. ) \;\; \mathrm{for} \;1 \le i < m,}{%
        r(i) = [ x_min + (i-1)*s,  x_min + i*s )   for 1<= i < m,}
  \deqn{r_m = \left[ \; x_{min} + (m-1)*s,\; x_{max} \; \right] ,}{%
        r(m) = [ x_min + (m-1)*s, x_max ] ,}

  where \eqn{x_{min}}{x_min} and \eqn{x_{max}}{x_max} are the minimum and the
  maximums of variable V respectively and
  \eqn{s = (x_{max} - x_{min}) / m}{s = (x_max - x_min) / m}.

}
\value{
  \item{tway.table}{two-way tables.}
  \item{total}{total number of data with corresponding code of variables.}
  \item{interval}{class interval for continuous and discrete explanatory
    variables.}
  \item{base.aic}{base_AIC.}
  \item{aic}{AIC's of single explanatory variables.}
  \item{aic.order}{list of explanatory variable numbers arranged in ascending
    order of AIC.}
  \item{nsub}{number of subsets of explanatory variables.}
  \item{subset}{list of subsets of explanatory variables in ascending order of
    AIC with the following components:

    \describe{
      \item{nv:}{number of explanatory variables,}
      \item{ncc:}{number of categories,}
      \item{aic:}{AIC's,}
      \item{exv:}{explanatory variables,}
      \item{vname:}{explanatory variable names.}
    }
  }
  \item{ctable}{contingency table constructed by the best subset and additional
    subsets if any variables is specified by \code{additional.output}.}
  \item{ctable.interval}{class interval for continuous and discrete explanatory
    variables in contingency table.}
  \item{caic}{AIC of subset of explanatory variables in contingency table.}
  \item{missing}{number of types of the missing values for each variable.}
}
\references{
  K.Katsura and Y.Sakamoto (1980) \emph{Computer Science Monograph, No.14,
  CATDAP, A Categorical Data Analysis Program Package.}
  The Institute of Statistical Mathematics.

  Y.Sakamoto (1985) \emph{Model Analysis of Categorical Data.} Kyoritsu Shuppan
  Co., Ltd., Tokyo. (in Japanese)

  Y.Sakamoto (1985) \emph{Categorical Data Analysis by AIC.} Kluwer Academic
  publishers.

  \emph{An AIC-based Tool for Data Visualization} (2015),
  \href{http://hdl.handle.net/10787/3614}{NTT DATA Mathematical Systems Inc.}
  (in Japanese)
}

\examples{
# Example 1 (medical data "HealthData")
# as additional output, contingency tables for explanatory variable sets
# c("aortic.wav","min.press") and c("ecg","age") are obtained.

data(HealthData)
catdap2(HealthData, c(2, 2, 2, 0, 0, 0, 0, 2), "symptoms",
        c(0., 0., 0., 1., 1., 1., 0.1, 0.), ,
        list(c("aortic.wav", "min.press"), c("ecg", "age")))

# Example 2 (Edgar Anderson's Iris Data)
# continuous response variable handling and the usage of Barplot2WayTable
# function to visualize the result in shape of stacked histogram.

data(iris)  
resvar <- "Petal.Width"
z <- catdap2(iris, c(0, 0, 0, -7, 2), resvar, c(0.1, 0.1, 0.1, 0.1, 0))
z

vname <- names(iris)
exvar <- c("Sepal.Length", "Petal.Length")
Barplot2WayTable(vname, resvar, exvar, z$tway.table, z$interval)

# Example 3  (in the case of a large number of variables)
data(HelloGoodbye)
pool <- rep(2, 56)

## using the default values of parameters pa1, pa2, pa3
## catdap2(HelloGoodbye, pool, "Isay", nvar = 10, print.level = 1, plot = 0) 
## Error : Working area for contingency table is too short, try pa1 = 12.

### According to the error message, set the parameter p1 at 12, then ..
catdap2(HelloGoodbye, pool, "Isay", nvar = 10, pa1 = 12, print.level = 1,
        plot = 0)

# Example 4 (HealthData with missing values)
data(MissingHealthData)
catdap2(MissingHealthData, c(2, 2, 2, 0, 0, 0, 0, 2), "symptoms",
        c(0., 0., 0., 1., 1., 1., 0.1, 0.), missingmark = 300)
}

\keyword{category}
