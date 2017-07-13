\name{GenAlg}
\alias{GenAlg}
\alias{newGeneration}
\alias{popDiversity}
\title{
  A generic Genetic Algorithm for feature selection
}
\description{
  These functions allow you to initialize (\code{GenAlg}) and iterate
  (\code{newGeneration}) a genetic algorithm to perform feature
  selection for binary class prediction in the context of gene
  expression microarrays or other high-throughput technologies.
}
\usage{
GenAlg(data, fitfun, mutfun, context, pm=0.001, pc=0.5, gen=1)
newGeneration(ga)
popDiversity(ga)
}
\arguments{
  \item{data}{
    The initial population of potential solutions, in the form of a data
    matrix with one individual per row.}
  \item{fitfun}{
    A function to compute the fitness of an individual solution. Must take
    two input arguments: a vector of indices into rows of the population
    matrix, and a \code{context} list within which any other items required
    by the function can be resolved. Must return a real number; higher values
    indicate better fitness, with the maximum fitness occurring at the optimal
    solution to the underlying numerical problem.}
  \item{mutfun}{
    A function to mutate individual alleles in the population. Must take two
    arguments: the starting allele and a \code{context} list as in the
    fitness function.}
  \item{context}{
    A list of additional data required to perform mutation or to compute
    fitness. This list is passed along as the second argument when
    \code{fitfun} and \code{mutfun} are called.}
  \item{pm}{
    A real value between \code{0} and \code{1}, representing the probability
    that an individual allele will be mutated.}
  \item{pc}{
    A real value between \code{0} and \code{1}, representing the probability
    that crossover will occur during reproduction.}
  \item{gen}{
    An integer identifying the current generation.}
  \item{ga}{
    An object of class \code{GenAlg}}
}
\value{
  Both the \code{GenAlg} generator and the \code{newGeneration} functions
  return a \code{\link{GenAlg-class}} object. The \code{popDiversity} function
  returns a real number representing the average diversity of the population.
  Here diversity is defined by the number of alleles (selected features) that
  differ in two individuals.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  P. Roebuck \email{proebuck@mdanderson.org}
}
\seealso{
  \code{\link{GenAlg-class}},
  \code{\link{GenAlg-tools}},
  \code{\link{maha}}.
}
\examples{
# generate some fake data
nFeatures <- 1000
nSamples <- 50
fakeData <- matrix(rnorm(nFeatures*nSamples), nrow=nFeatures, ncol=nSamples)
fakeGroups <- sample(c(0,1), nSamples, replace=TRUE)
myContext <- list(dataset=fakeData, gps=fakeGroups)

# initialize population
n.individuals <- 200
n.features <- 9
y <- matrix(0, n.individuals, n.features)
for (i in 1:n.individuals) {
  y[i,] <- sample(1:nrow(fakeData), n.features)
}

# set up the genetic algorithm
my.ga <- GenAlg(y, selectionFitness, selectionMutate, myContext, 0.001, 0.75)

# advance one generation
my.ga <- newGeneration(my.ga)

}
\keyword{optimize}

