\name{MWLasso-package}
\alias{MWLasso-package}
\alias{MWLasso}
\docType{package}
\title{
A penalized Moving-Window Lasso (MWLasso) method for genome-wide association studies.
}
\Version: \tab 1.0\cr
\Date: \tab 2016-07-25\cr

\author{
Minli Bao <minli-bao@uiowa.edu>
}

\tabular{lll}{
Package: \tab MWLasso\cr
Type: \tab Package\cr
License: \tab GPL\cr
}


\description{
Find coefficients for a regression model with moving-window Lasso penalty (MWLasso). The model is marginalized so that missing values can be accommendated.
}

\usage{
MWLasso(X, Y, lambda, eta, d, method='linear', epson = 1e-10, M = 100)
MW_parameters(X, Y, d, number, gamma2, method="linear", epson = 1e-10, M = 100)
}

\arguments{
\item{X}{ The design matrix which can include missing values.}
\item{Y}{The response variable.}
\item{lambda}{The L1 penalty tuning parameter.}
\item{eta}{The moving-window penalty tuning parameter.}
\item{d}{The size of the moving-window penalty.}
\item{method}{The regression method. Method can be linear or logistic.}
\item{number}{The pre-determined number of non-zero coefficients.}
\item{gamma2}{The proportion of the tuning parameter for moving-window part and L1 penalty part. gamma2 is the ratio of lambda/(lambda+eta).}
\item{epson}{Convergence criterion. The iteration will stop if the relative change is smaller than epson.}
\item{M}{The maximum number of iterations.}
}
\details{
The function minimizes 1/(2n)*MLS + lambda*L1 + eta/(2(d-1))*MW. Here MLS is the marginalized least squares, L1 is the L1 penalty in Lasso, and MW is the moving-window penalty.
}
\value{
MWLasso returns:
\item{beta}{The coefficients estimates.}
MW_parameters returns:
\item{lambda}{The tuning parameter for L1 penalty.}
\item{eta}{The tuning parameter for moving-window penalty.}
Thus under such lambda and eta, MWLasso will return beta which satisfies the pre-determined number of non-zero coefficients.
}
\references{
Genome-wide association studies using a penalized moving-window regression. Minli Bao and Kai Wang. Submitted for review.
}

\keyword{ package }

\examples{
n = 100
p = 200
X = matrix(rnorm(n*p, mean = 0, sd = 1), nrow = n, ncol = p)
X[2,1] = NA; X[80:90,1]= NA; 
Y =  rnorm(n, mean = 0, sd = 1)

number = 20
gamma2 = 0.1
d = 4
params = MW_parameters(X, Y, d, number, gamma2, method="linear", epson = 1e-10, M = 100)
lambda = params[1]
eta = params[2]
beta.hat = MWLasso(X, Y, lambda, eta, d, method='linear')
}
