\name{maha}
\Rdversion{1.1}
\alias{maha}
\title{
Compute the (squared) Mahalanobis distance between two groups of vectors
}
\description{
The Mahalanobis distance between two groups of vectors 
}
\usage{
maha(data, groups, method = "mve")
}

\arguments{
  \item{data}{A matrix with columns representing features (or variables)
    and rows representing independent samples}
  \item{groups}{
    A factor or logical vector with length equal to the number of rows
    (samples) in the \code{data} matrix}
  \item{method}{
    A character string determining the method that should be used to
    estimate the covariance matrix. The default value of "\code{mve}"
    uses the \link[MASS]{cov.mve} function from the MASS package. The
    other valid option is "var", which uses the \code{\link{var}}
    function from the standard \code{stats} package.}
}
\details{
  The Mahalanobis distance between two groups of vectors is the distance
  between their centers, computed in the equivalent of a principal
  component space that accounts for different variances.
}
\value{
  Returns a numeric vector of length 1.
}
\references{
  Mardia, K. V. and Kent, J. T. and Bibby, J. M.\cr
  \emph{Multivariate Analysis}.\cr
  Academic Press, Reading, MA {1979}, pp. 213--254.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  P. Roebuck \email{proebuck@mdanderson.org}
}
\seealso{
  \code{\link[MASS]{cov.mve}},
  \code{\link[stats]{var}}
}
\examples{
nFeatures <- 40
nSamples <- 2*10
dataset <- matrix(rnorm(nSamples*nFeatures), ncol=nSamples)
groups <- factor(rep(c("A", "B"), each=10))
maha(dataset, groups)
}
\keyword{multivariate}

