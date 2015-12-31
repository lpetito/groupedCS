\name{PAVA for Current Status Data}
\alias{pava.cs}
\alias{pava.cs.mc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variations on PAVA
}
\description{
These functions are the variations on the pool-adjacent-violators algorithm needed to compute the NPMLE for current status data (\code{pava.cs}) and misclassified current status data (\code{pava.cs.mc}).
}
\usage{
pava.cs(Cs, initial)
pava.cs(Cs, initial, alpha, beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Cs}{
Observation times for each individual
}
  \item{initial}{
The test result for each individual (true if \code{alpha} and \code{beta} are 1, misclassified otherwise)
}
  \item{alpha}{
Sensitivity. Default is 1.
}
  \item{beta}{
Specificity. Default is 1.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The output of the PAV (or adjusted PAV) algorithm. This is just a list of the predicted cumulative probabilities at each value of C. 
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data <- gen.data.weibull.unif(100, 2, 4, 25, 25, 0.9, 0.9)
with(data, pava.cs(data$Cs, data$delta.ind))
with(data, pava.cs.mc(data$Cs, data$y.ind, alpha=0.9, beta=0.9))
}
\keyword{survival}
