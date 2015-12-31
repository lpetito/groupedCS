\name{Simulation function - fixed censoring}
\alias{simulation.fixed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation function - fixed censoring
}
\description{
This function does one replication of the simulation for the supplemental materials section of the paper from data generated with fixed censoring times and a user-specified true event probability. It returns a description of the misclassification of both the individual and group tests, the results from the appropriate PAVA, the results from the hybrid EM-PAV algorithm for grouped tests, and the number of iterations the hybrid EM-PAV algorithm takes to converge
}
\usage{
simulation.fixed(n, k, Cs, true.F, alpha, beta, t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
number of individuals
}
  \item{k}{
grouping size
}
  \item{Cs}{
a vector of the observed censoring times
}
  \item{true.F}{
a vector of event probabilities at each one of the \code{Cs}. This vector must be the same length as \code{Cs}.
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
The same list as returned by \code{simulation.random}
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
\code{\link{hybrid.em.pav}}, \code{\link{pava.cs.mc}}, \code{\link{gen.data.fixed}}
}
\examples{
simulation(100, 2, 1:10, seq(0.05, 0.5, 0.05), 0.95, 0.95, 0.01)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
