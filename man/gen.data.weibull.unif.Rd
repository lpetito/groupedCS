\name{Generate data with Weibull events and Uniform Censoring}
\alias{gen.data.weibull.unif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data Generation
}
\description{
This function generates grouped, misclassified current status data with event times from a Weibull distribution with parameters shape and scale (as defined in the stats package) and observation times from a Uniform distribution that allows roughly quantile percent of the events to be observed.
}
\usage{
gen.data.weibull.unif(n, k, shape, scale, quantile = 0.99, alpha=1, beta=1)
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
  \item{alpha}{
Sensitivity: probability of a positive test results given that the individual is truly diseased (or that the group contains at least one person who is truly diseased). Default is 1 - no misclassification
}
  \item{beta}{
Specificity: probability of a negative test results given that the individual is truly not diseased (or that the group contains noone who is truly diseased). Default is 1 - no misclassification
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This returns a data frame with the following columns:
1) groups: group identifier
2) Ts: true individual event times
3) Cs: individual observation times 
4) delta.ind: indicator of event (1) or censoring (0) - true test result
5) y.ind: misclassified test result 
6) delta.group: true group test result, indicator that at least one individual had delta.ind==1
7) y.group: misclassified group test result
8) initial.p: initial values for the EM-PAV hybrid algorithm
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Lucia Petito
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gen.data.fixed}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Generate data for 10 people with group sizes of 2 from a Weibull(4, 25) distribution where at most 50\% of the population can be diseased, from a test with a 10\% false positive and negative rate. 
data <- gen.data.weibull.unif(10, 2, 4, 25, 0.5, 0.9, 0.9)
data

#Redo the above scenario, but use 1,000 people and summarize censoring and misclassification rates.
data <- gen.data.weibull.unif(1000, 2, 4, 25, 0.5, 0.9, 0.9)
head(data)
#Number of individual events
with(data, xtabs(~delta.ind))
#Number of group events
with(data, xtabs(~delta.group))
#Summary of misclassification of individual events
with(data, xtabs(~delta.ind + y.ind))
#Summary of misclassification of group events
with(data, xtabs(~delta.group + y.group))

}

\keyword{survival}
\keyword{datagen}
