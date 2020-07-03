

gibbsparam <- function(data, hyper, var, sigm, nbeps, qq) {
  # gibbs step for the parmeter, prior induced by the hyperparameter
  p <- length(data)
  thetc <- numeric(p)
  dists <- 0
  for (i in 1:p) {
    thettest <- rnorm(nbeps, hyper, var)
    test <- rowMeans(matrix(rnorm(qq * nbeps, thettest, sigm), nrow = nbeps))
    dist <- abs(test - data[i])
    thetc[i] <- thettest[which.min(dist)]
    dists <- dists + min(dist)
  }
  return(list(thetc, dists))
}



gibbshyper <- function(thet, nbeps2, var) {
  # gibbs step for the hyperparameter, with uniform prior
  res <- runif(nbeps2, -4, 4)
  test <- rowMeans(matrix(rnorm(length(thet) * nbeps2, res, var), ncol = length(thet)))
  dist <- abs(test - mean(thet))
  return(list(res[which.min(dist)], min(dist)))
}

gibbstot <- function(data, thetini, hyperini, sigm, var, nbeps1, nbeps2, nbpts, qq) {
  # full function
  reshyper <- rep(NA, nbpts + 1)
  resparam <- matrix(NA, ncol = nbpts + 1, nrow = length(thetini))
  reshyper[1] <- hyperini
  resparam[, 1] <- thetini
  resdist <- rep(NA, nbpts)
  for (i in 2:(nbpts + 1)) {
    resdist[i - 1] <- 0
    VV <- gibbsparam(data, reshyper[i - 1], var, sigm, nbeps1, qq)
    resparam[, i] <- VV[[1]]
    resdist[i - 1] <- resdist[i - 1] + VV[[2]]
    WW <- gibbshyper(resparam[, i], nbeps2, var)
    reshyper[i] <- WW[[1]]
  }
  return(list(reshyper, resparam, resdist))
}


# for comparison, simple ABC
# here, ntot=nbpts*(qq*nbeps1 + nbeps2)
ABCsimple <- function(data, vari, sigma, nbpts, ntot, qq) {
  reshyper <- rep(NA, ntot)
  resparam <- matrix(NA, ncol = ntot, nrow = length(data))
  dist <- matrix(NA, nrow = length(data), ncol = ntot)
  for (i in 1:ntot) {
    reshyper[i] <- runif(1, -4, 4)
    resparam[, i] <- rnorm(length(data), reshyper[i], vari)
    dist[i] <- sum(abs(data - rowMeans(
      matrix(rnorm(length(data) * qq, rep(resparam[, i], qq), sigma), nrow = length(data))
    )))
  }
  # compute the distance in the end
  v <- order(dist)[1:nbpts]
  return(list(reshyper[v], resparam[, v], dist[v[nbpts]]))
}


# estimation of the variance of the estimator
resgibbs <- function(data, thetini, hyperini, sigm, var, nbeps1, nbeps2, nbpts, nbvar, qq) {
  res1 <- rep(NA, nbvar)
  res2 <- matrix(NA, ncol = nbvar, nrow = length(thetini))
  for (i in 1:nbvar)
  {
    U <- gibbstot(data, thetini, hyperini, sigm, var, nbeps1, nbeps2, nbpts, qq)
    res1[i] <- mean(U[[1]][2:(nbpts + 1)])
    res2[, i] <- rowSums(U[[2]][, 2:(nbpts + 1)]) / nbpts
  }
  return(list(var(res1), mean(res1), apply(res2, 1, var), apply(res2, 1, mean)))
}

ressimple <- function(data, sigm, var, nbpts, ntot, nbvar, qq) {
  res1 <- rep(NA, nbvar)
  res2 <- matrix(NA, ncol = nbvar, nrow = length(data))
  for (i in 1:nbvar)
  {
    U <- ABCsimple(data, var, sigm, nbpts, ntot, qq)
    res1[i] <- mean(U[[1]])
    res2[, i] <- rowSums(U[[2]]) / nbpts
  }
  return(list(var(res1), mean(res1), apply(res2, 1, var), apply(res2, 1, mean)))
}

# Gibbs exact
gibbsexact <- function(data, k, sigm, var, qq) {
  R1 <- matrix(NA, ncol = 20, nrow = k)
  R2 <- rep(NA, k)
  R2[1] <- rnorm(1, -4, 4)
  R1[1, ] <- rnorm(20, (qq * data / var + R2[1] / sigm) / (1 / sigm + qq / var), 1 / (1 / sigm + qq / var))
  for (i in 2:k)
  {
    cand <- runif(1, -4, 4)
    u <- runif(1, 0, 1)
    if (u < (prod(dnorm(R1[i - 1, ], cand, var))) / (prod(dnorm(R1[i - 1, ], R2[i - 1], var)))) {
      R2[i] <- cand
    } else {
      R2[i] <- R2[i - 1]
    }
    R1[i, ] <- rnorm(20, (qq * data / var + R2[i] / sigm) / (1 / sigm + qq / var), 1 / (1 / sigm + qq / var))
  }
  return(list(R1, R2))
}


# un SMC maison


SMCmaison <- function(npart, M, itermax, epstarget, alph, data, prior, model) {
  par <- sapply(1:npart, function(x) {
    prior$simu()
  })
  stats <- sapply(1:npart, function(x) {
    sapply(1:M, function(y) {
      model(par[, x], data)
    })
  })
  j <- 1
  pds <- rep(1 / npart, npart)
  ESS <- npart
  eps <- max(stats)
  histeps <- eps
  histpart <- matrix(ncol = npart, nrow = 2 * itermax)
  histpds <- matrix(ncol = npart, nrow = itermax)
  histpds[1, ] <- pds
  cb <- rep(0, itermax)
  while (j < itermax && eps > epstarget) {
    VV <- changementeps(stats, eps, pds, alph)
    eps <- VV[[1]]
    pds <- VV[[2]]
    ESS <- VV[[3]]
    histeps <- c(histeps, eps)
    histpart[2 * (j - 1) + 1, ] <- par[1, ]
    if (ESS < npart / 2) {
      sd <- 2 * cov(t(par[, which(pds != 0)]))
      if (any(diag(sd) == 0)) {
        sd <- sd + diag(.1, nrow(par))
      }
      part <- par
      statst <- stats
      for (i in 1:npart) {
        NON <- TRUE
        kk <- 1
        while (NON & kk < 10000) {
          qui <- sample(1:npart, 1, prob = pds, rep = TRUE)
          para <- par[, qui]
          statsa <- stats[, qui]
          part <- para + rmvnorm(1, sigma = sd)
          statst <- sapply(1:length(statsa), function(x) {
            model(part, data)
          })
          u <- runif(1)
          kk <- kk + 1
          if (u < (sum(statst <= eps) / sum(statsa <= eps) * prior$density(part) / prior$density(para))) {
            par[, i] <- part
            stats[, i] <- statst
            NON <- FALSE
          }
        }
        if (kk == 10000) {
          qui <- sample(1:npart, 1, prob = pds, rep = TRUE)
          para <- par[, qui]
          statsa <- stats[, qui]
        }
        cb[j] <- cb[j] + kk
      }
      pds <- rep(1 / npart, npart)
      ESS <- npart
      histpart[2 * j, ] <- par[1, ]
    }
    VV <- pasnoyau(par, eps, stats, pds, data, prior, model)
    par <- VV[[1]]
    histpds[1, ] <- pds
    stats <- VV[[2]]
    print(paste(j, eps))
    j <- j + 1
  }
  # on finit par un resampling pour avoir quelque chose de propre
  qui <- sample(1:npart, prob = pds, rep = TRUE)
  par <- par[, qui]
  pds <- rep(1 / npart, npart)
  ESS <- npart
  stats <- stats[, qui]
  return(list(par, pds, stats, histeps, histpart, histpds, cb))
}


SMCmaisonold <- function(npart, M, itermax, epstarget, alph, data, prior, model) {
  par <- sapply(1:npart, function(x) {
    prior$simu()
  })
  stats <- sapply(1:npart, function(x) {
    sapply(1:M, function(y) {
      model(par[, x], data)
    })
  })
  j <- 1
  pds <- rep(1 / npart, npart)
  ESS <- npart
  eps <- max(stats)
  histeps <- eps
  while (j < itermax && eps > epstarget) {
    VV <- changementeps(stats, eps, pds, alph)
    eps <- VV[[1]]
    pds <- VV[[2]]
    ESS <- VV[[3]]
    histeps <- c(histeps, eps)
    if (ESS < npart / 2) {
      qui <- sample(1:npart, prob = pds, rep = TRUE)
      par <- par[, qui]
      pds <- rep(1 / npart, npart)
      ESS <- npart
      stats <- stats[, qui]
    }
    VV <- pasnoyau(par, eps, stats, pds, data, prior, model)
    par <- VV[[1]]
    stats <- VV[[2]]
    print(paste(j, eps))
    j <- j + 1
  }
  # on finit par un resampling pour avoir quelque chose de propre
  qui <- sample(1:npart, prob = pds, rep = TRUE)
  par <- par[, qui]
  pds <- rep(1 / npart, npart)
  ESS <- npart
  stats <- stats[, qui]
  return(list(par, pds, stats, histeps))
}



# par : matrice des paramètres des particules actuelles, en colonne les coordonnées des paramètres
# simu : matrice des distances entre obs et simus
pasnoyau <- function(par, eps, stats, pds, data, prior, model) {
  sd <- 2 * cov(t(par))
  if (any(diag(sd) == 0)) {
    sd <- sd + diag(.1, nrow(par))
  }
  part <- par
  statst <- stats
  for (i in 1:ncol(par)) {
    if (pds[i] > 0) {
      VV <- chgtpar(par[, i], stats[, i], sd, eps, data, prior, model)
      part[, i] <- VV[[1]]
      statst[, i] <- VV[[2]]
    }
  }
  return(list(part, statst))
}

changementeps <- function(stats, eps, pds, alph) {
  posseps <- sort(stats[stats <= eps], decreasing = TRUE)
  pdst <- pds
  epst <- eps
  test <- FALSE
  k <- 1
  newpds <- pds
  ESS <- 1 / (sum(pds^2))
  ESSt <- ESS
  while (k < (length(posseps) - 1) && !test) {
    neweps <- posseps[k]
    for (i in 1:length(pds)) {
      if (newpds[i] != 0) {
        newpds[i] <- pds[i] * sum(stats[, i] <= neweps) / sum(stats[, i] <= eps)
      }
    }
    newpds <- newpds / sum(newpds)
    newESS <- 1 / sum(newpds^2)
    if (newESS < alph * ESS) {
      test <- TRUE
      pdst <- newpds
      epst <- neweps
      ESSt <- newESS
    }
    k <- k + 1
  }
  return(list(epst, pdst, ESSt))
}


chgtpar <- function(par, stats, sd, eps, data, prior, model) {
  part <- par + rmvnorm(1, sigma = sd)
  statst <- sapply(1:length(stats), function(x) {
    model(part, data)
  })
  u <- runif(1)
  if (u < (sum(statst <= eps) / sum(stats <= eps) * prior$density(part) / prior$density(par))) {
    return(list(part, statst))
  } else {
    return(list(par, stats))
  }
}
