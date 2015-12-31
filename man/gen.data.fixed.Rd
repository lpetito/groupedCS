\name{Generate data from a fixed distribution on a fixed finite number of censoring times}
\alias{gen.data.fixed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fixed data generation
}
\description{
This function generates grouped, misclassified current status data with a finite number of fixed censoring times. Probability of an event is determined from a user-provided "distribution" - only requirement is the vector \code{true.F} must be monotonically increasing and between 0 and 1.
}
\usage{
gen.data.fixed.Rd(n, k, Cs, true.F, alpha=1, beta=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{number of individuals}
  \item{k}{grouping size}
  \item{Cs}{a vector of the observed censoring times}
  \item{true.F}{a vector of event probabilities at each one of the \code{Cs}. This vector must be the same length as \code{Cs}.}
  \item{alpha}{Sensitivity: probability of a positive test results given that the individual is truly diseased (or that the group contains at least one person who is truly diseased). Default is 1 - no misclassification}
  \item{beta}{Specificity: probability of a negative test results given that the individual is truly not diseased (or that the group contains noone who is truly diseased). Default is 1 - no misclassification}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function returns a data frame with the following columns:
\tabular{ll}{
Cs: \tab individual observation times \cr
groups: \tab group identifier\cr
initial.p: \tab initial values for the EM-PAV hybrid algorithm\cr
delta.ind: \tab indicator of event (1) or censoring (0) - true test result\cr
y.ind: \tab misclassified test result\cr
delta.group: \tab true group test result, indicator that at least one individual had delta.ind==1\cr
y.group: \tab misclassified group test result\cr
}
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
\code{\link{gen.data.weibull.unif}}
}
\examples{
#Generate data on 10 individuals with 5\% misclassification rates each
data <- gen.data.fixed(10, 2, 1:5, seq(0.1, 0.5, 0.1), 0.95, 0.95)
data

#Now examine generated data in 1,000 individuals 
data <- gen.data.fixed(1000, 2, 1:5, seq(0.1, 0.5, 0.1), 0.9, 0.9)
#Look at true individual test results
with(data, xtabs(~Cs + delta.ind))
#Look at misclassification
with(data, xtabs(~delta.ind + y.ind + Cs))

#Do the same in the grouped tests 
with(data, xtabs(~Cs + delta.group))
with(data, xtabs(~delta.group + y.group + Cs))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
