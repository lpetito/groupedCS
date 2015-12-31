\name{Simulation function}
\alias{simulation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation function
}
\description{
This function does one replication of the simulation for the paper. It returns a description of the misclassification of both the individual and group tests, the results from the appropriate PAVA, the results from the hybrid EM-PAV algorithm for grouped tests, and the number of iterations the hybrid EM-PAV algorithm takes to converge
}
\usage{
simulation(n, k, shape, scale, quantile, x, alpha=1, beta=1, t=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
number of individuals
}
  \item{k}{
grouping size
}
  \item{shape}{
shape for the Weibull distribution
}
  \item{scale}{
scale for the Weibull distribution (defaults to 1)
}
  \item{quantile}{
The maximum probability of the event in the population (default is 0.99)
}
  \item{x}{
The maximum value for censoring. This should correspond with the quantile argument.
}
  \item{alpha}{
Sensitivity: probability of a positive test results given that the individual is truly diseased (or that the group contains at least one person who is truly diseased). Default is 1 - no misclassification
}
  \item{beta}{
Specificity: probability of a negative test results given that the individual is truly not diseased (or that the group contains noone who is truly diseased). Default is 1 - no misclassification
}
  \item{t}{
threshold for convergence (default is 0.01)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{desc.ind}{Table with description of the misclassification of the individual tests}
\item{desc.group}{Table with description of the misclassification of the group tests}
\item{num.it}{Number of iterations for the hybrid EM-PAV algorithm to converge}
\item{ind.result}{Result from appropriate PAV algorithm (\code{pava.cs} if \code{alpha} = \code{beta} = 1, \code{pava.cs.mc} otherwise)}
\item{group.result}{Result from hybrid EM-PAV algorithm, see function \code{hybrid.em.pav} for details}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Lucia Petito
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{hybrid.em.pav}}, \code{\link{pava.cs.mc}}, \code{\link{gen.data.weibull.unif}}
}

\examples{

simulation(100, 2, 4, 25, 0.25, 14, 0.95, 0.95, 0.01)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

