\name{catdap2}
\alias{catdap2}

\title{Categorical Data Analysis Program Package 02}
\description{
  Search for the best single explanatory variable and detect the best subset of explanatory variables.
}
\usage{
catdap2(data, pool=NULL, response.name, accuracy=NULL, nvar=NULL,
        explanatory.names=NULL, pa1=1, pa2=4, pa3=10, print.level=0, plot=1)
}
\arguments{
  \item{data}{data matrix with variable names on the first row.}
  \item{pool}{the ways of pooling to categorize each variable must be specified by integer parameters.
    \tabular{rl}{
      (-m) < 0 : \tab m-bin histogram is employed to describe the distribution of continuous response variable. \cr
      \tab (this option is valid only for the response variable) \cr 
      = 0 : \tab equally spaced pooling via a top-down algorithm \cr
      = 1 : \tab unequally spaced pooling via a bottom-up algorithm (default) \cr
      = 2 : \tab no pooling for discrete variables \cr
    }}
  \item{response.name}{variable name of the response variable.}
  \item{accuracy}{minimum width for the discretization for each variable.}
  \item{nvar}{number of variables to be retained for the analysis of multidimensional tables. Default is the number of variables in \code{data}.}
  \item{explanatory.names}{set of variable names  for additional analyses.}
  \item{pa1, pa2, pa3}{control parameter for size of the working area. If error message is output, please change the value of parameter according to it.}
  \item{print.level}{this argument determines the level of output printing. The default value of \code{0} means that lists of "AIC's of the models with k explanatory variables (k=1,2,...)" are printed. A value of \code{1} means that those lists are not printed and "Summary of subsets of explanatory variables" within the top 30 is listed.}
  \item{plot}{split directions for each level of the mosaic.
    \tabular{rl}{
      = 0 : \tab no plot \cr
      = 1 : \tab horizontal (default) \cr
      = 2 : \tab alternating directions, beginning with a vertical split \cr
    }}
}
\details{
  This function is a R-function style clone of Sakamoto's CATDAP-02 program
  for categorical data analysis. CATDAP-02 can be used to search for the best subset
  of explanatory variables which have the most effective information on a specified
  response variable. Continuous explanatory variables could be explanatory variables.
  In that case CATDAP-02 searches for optimal categorization of continuous values.

  The basic statistic adopted is obtained by the application of the
  statistic AIC to the models.

  \eqn{E} denotes the response variable and \eqn{F} denotes candidate explanatory 
  variable, and their cell frequencies by \eqn{n_E(i) (i\in(E))}{n_E(i) (i in (E))}
  and \eqn{n_F(j) (j\in(F))}{n_F(j) (j in (F))}. The cross frequency is denoted by
  \eqn{n_F,F(i,j)}{n_F,F(i,j)} \eqn{(i,j\in(E,F))}{(i,j in (E,F))}. 
  To measure the strength of dependence of a specific set of response
  variables E on the explanatory variable F, we use the
  following statistic:
    \deqn{AIC(E;F)=-2\sum_{i,j\in(E,F)}n_E,F(i,j)\ ln\{n_E,F(i,j)/(n_F(j)\} + 2(C_E-1)C_F,\ \ (1)}{AIC(E;F)=-2 \sum _{i,j in (E,F)}n_E,F(i,j) ln{n_E,F(i,j)/(n_F(j)} + 2(C_E-1)C_F,  (1)}

  where \eqn{C_E} and \eqn{C_F} denote the total number of categories
  of the corresponding sets of variables, respectively. 

  The selection of the best subset of explanatory variables is
  realized by the search for \eqn{F} which gives the minimum \eqn{AIC(E;F)}.

  In case of \eqn{F=\phi}, the formula (1) reduces to

    \deqn{AIC(E;\phi)=-2\sum_{i\in(E)}n_E(i)\ ln\{n_E(i)/n\}+ 2(C_E-1),\ \ \ (2)}{AIC(E;\phi)=-2 \sum _{i in (E)}n_E(i) ln{n_E(i)/n}+ 2(C_E-1),   (2)}

  Here it is assumed that \eqn{C_\phi=1} and \eqn{n_\phi(1)=n}.

  Sakamoto's original CATDAP outputs \eqn{AIC(E;F) - AIC(E;\phi)} 
  as the AIC value instead of \eqn{AIC(E;F)}. 
  By this way the positive value of AIC indicates that the variable 
  \eqn{F} is judged to be useless as the explanatory variable of 
  the \eqn{E}.

  On the other hand, this policy make impossible to compare the
  goodness of the CATDAP model with other models, logit models 
  for example. 

  Considering the convenience of users, present "R version CATDAP"
  provides not only \eqn{AIC = AIC(E;F) - AIC(E;\phi)}, but 
  \eqn{AIC(E;\phi)}, either. The latter value is given as base_AIC in the output.

  Users could recover \eqn{AIC(E;F)} by adding AIC and base_AIC.
}
\value{
  \item{tway.table}{two-way tables.}
  \item{total}{total number of data with corresponding code of variables.}
  \item{interval}{class interval for continuous and discrete explanatory variables.}
  \item{base.aic}{base_AIC.}
  \item{aic}{AIC's of explanatory variables.}
  \item{aic.order}{list of explanatory variable numbers arranged in ascending order of AIC.}
  \item{nsub}{number of subsets of explanatory variables.}
  \item{nv}{number of explanatory variables in a subset.}
  \item{ncc}{number of categories in a subset of explanatory variables.}
  \item{aaic}{AIC's of subsets of explanatory variables.}
  \item{cexp}{list of subsets of explanatory variables.} 
  \item{ctable}{contingency table constructed by the best subset and additional subsets if any variables is specified by "explanatory.names".}
  \item{ctable.interval}{class interval for continuous and discrete explanatory variables in contingency table.}
  \item{caic}{AIC of subset of explanatory variables in contingency table.}
}
\references{
  K.Katsura and Y.Sakamoto (1980) \emph{Computer Science Monograph, No.14, CATDAP, A Categorical Data Analysis Program Package.} The Institute of Statistical Mathematics.

  Y.Sakamoto (1985) \emph{Model Analysis of Categorical Data.} Kyoritsu Shuppan Co., Ltd., Tokyo. (in Japanese)

  Y.Sakamoto (1985) \emph{Categorical Data Analysis by AIC.} Kluwer Academic publishers.

  \emph{An AIC-based Tool for Data Visualization} (2015), \href{http://ismrepo.ism.ac.jp/dspace/browse?type=author&value=NTT+DATA+Mathematical+Systems+Inc.}{NTT DATA Mathematical Systems Inc.} (in Japanese)
}

\examples{
# Example 1 (Application to medical data "HealthData")
# as additional analysis, a contingency table for explanatory variables "aortic.wav"
# and "min.press" is obtained.
data(HealthData)
catdap2(HealthData, c(2, 2, 2, 0, 0, 0, 0, 2), "symptoms",
        c(0., 0., 0., 1., 1., 1., 0.1, 0.), , c("aortic.wav", "min.press"))

# Example 2 (continuous response variable is detected)
# Edgar Anderson's Iris Data
# "Sepal.Length" "Sepal.Width"  "Petal.Length" "Petal.Width"  "Species" 
data(iris)  
catdap2(iris, c(0,0,0,-7,2), "Petal.Width", c(0.1, 0.1, 0.1, 0.1, 0))

# Example 3  (in the case of a large number of variables)
data(HelloGoodbye)
pool <- rep(2,56)

## using the default values of parameters pa1,pa2,pa3
## catdap2(HelloGoodbye, pool, "Isay", nvar=10, print.level=1, plot=0) 
## error : Working area for contingency table is too short, try pa1=12.

### According to the error message, set the parameter p1 at 12, then ..
catdap2(HelloGoodbye, pool, "Isay", nvar=10, pa1=12, print.level=1, plot=0)
}

\keyword{category}
