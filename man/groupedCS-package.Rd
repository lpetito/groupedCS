\name{groupedCS-package}
\alias{groupedCS-package}
\alias{groupedCS}
\docType{package}
\title{
groupedCS
}
\description{
This package contains all of the code needed to execute all of the simulations presented in the working paper "Misclassified Group Tested Current Status Data". It includes data generating functions, the hybrid EM-PAV algorithm, and the appropriate version of the pool-adjacent-violators algorithm for perfectly classified and misclassified current status data.
}
\details{
\tabular{ll}{
Package: \tab groupedCS\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-12-31\cr
License: \tab GPL-3\cr
}
}
\author{
Author and Maintainer: Lucia Petito <lucia.petito@gmail.com>
}
\seealso{
\code{\link{hybrid.em.pav}}, \code{\link{pava.cs.mc}}, \code{\link{gen.data.weibull.unif}}
}
\examples{
#For a detailed example of how the code is used, see vingette.

#To run one iteration of the code for uniform censoring times, use
sim.r <- simulation.random(100, 2, 4, 25, 0.25, 14, 0.95, 0.95, 0.01)
#See output
sim.r

#Do the same for fixed censoring times
sim.f <- simulation.fixed(100, 2, 1:5, seq(0.1, 0.5, 0.1), 0.95, 0.95, 0.01)
#See output
sim.f
}
\keyword{ package }