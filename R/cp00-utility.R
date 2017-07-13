# Copyright (C) Kevin R. Coombes, 2007-2012

# This routine computes the Mahalanobis distance between the centers of
# two groups. We sometimes use this as the fitness function in our genetic
# algorithm.
maha <- function(data, groups, method='mve') {
  if (is.logical(groups)) {
    og <- groups
    groups <- rep('B', length(groups))
    groups[og] <- 'A'
  }
  if (!is.factor(groups)) {
    groups <- factor(groups)
  }
  xa <- apply(data[groups==levels(groups)[[1]],], 2, mean)
  xb <- apply(data[groups==levels(groups)[[2]],], 2, mean)
  if (method == 'mve') {
    covar <- cov.mve(data)$cov
  } else {
    covar <- var(data)
  }
  t(xa - xb) %*% solve(covar) %*% (xa - xb)
}

