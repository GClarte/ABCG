library(ggplot2)
library(ggpubr)
library(smfsb)
library(EasyABC)
library(mvtnorm)

source("functions_hierarchical.R")

# generation of the dataset
a <- 2
# on code les variances en dur pour éviter les problèmes
vari <- 1
sigma <- 1
Nk <- 10
thetreel <- c(2.774149, 2.438900, 3.538051, 3.134760, 3.192536, 2.134164, 2.460961, 2.785255, 1.821516, 1.999794, 2.147105, 2.627481, 2.273351, 1.829135, 1.095750, 2.192478, 1.951274, 2.781893, 1.795393, 2.563236)
obs <- c(2.5666002, 2.5808230, 3.4319618, 3.4865059, 3.5939341, 1.7859019, 2.3653864, 2.7560831, 2.2540691, 2.0810668, 2.3309626, 2.1923769, 2.8244147, 2.3624706, 0.8032528, 1.7128994, 2.1605168, 3.5554632, 1.4499120, 2.3010078)

# pour SMCmaison

rpriortot <- function() {
  alpha <- runif(1, -4, 4)
  thet <- rnorm(20, alpha)
  return(c(alpha, thet))
}

dpriortot <- function(x) {
  if (abs(x[1]) < 4) {
    return(1 / 4 * prod(dnorm(x[2:21], x[1])))
  } else {
    return(0)
  }
}

monprior <- list(density = dpriortot, simu = rpriortot)

model <- function(pars, data) {
  return(sum(abs(rowMeans(matrix(rnorm(20 * 10, rep(pars[2:21], 10), 1), nrow = 20)) - data)))
}

# pour plusieurs valeurs des paramètres
# le nombre de simulations du modèle vaut 100*1000
resG100 <- gibbstot(obs, rep(0, 20), 0, sigma, vari, 100, 100, 1000, 10)
resS100 <- ABCsimple(obs, vari, sigma, 1000, 100 * 1000, 10)

resG10 <- gibbstot(obs, rep(0, 20), 0, sigma, vari, 10, 10, 1000, 10)

resG1 <- gibbstot(obs, rep(0, 20), 0, sigma, vari, 1, 1, 1000, 10)
resG5 <- gibbstot(obs, rep(0, 20), 0, sigma, vari, 5, 5, 1000, 10)
resG10 <- gibbstot(obs, rep(0, 20), 0, sigma, vari, 10, 10, 1000, 10)
resG1000 <- gibbstot(obs, rep(0, 20), 0, sigma, vari, 1000, 1000, 1000, 10)

exact <- gibbsexact(obs, 10000, sigm, var, 10)


# pour notre SMC maison

resP100 <- SMCmaison(1000, 20, 1000, resG100[[3]], .9, obs, monprior, model)


# pour superposition à coût constant
hyper <- numeric()
par1 <- numeric()
par2 <- numeric()
par3 <- numeric()
nbr <- numeric()
type <- numeric()
exact <- gibbsexact(obs, 1000, sigma, vari, 10)


N <- 100
for (i in 1:N) {
  A <- gibbstot(obs, rep(0, 20), 0, sigma, vari, 50, 50, 50, 10)
  B <- ABCsimple(obs, vari, sigma, 50, 2500, 10)
  C <- SMCmaison(50, 2, 25, 0, .9, obs, monprior, model)
  hyper <- c(hyper, A[[1]][5:50], B[[1]], C[[1]][1, ])
  par1 <- c(par1, A[[2]][1, 5:50], B[[2]][1, ], C[[1]][2, ])
  par2 <- c(par2, A[[2]][2, 5:50], B[[2]][2, ], C[[1]][3, ])
  par3 <- c(par3, A[[2]][3, 5:50], B[[2]][3, ], C[[1]][4, ])
  nbr <- c(nbr, rep(i, length(A[[1]][5:50]) + length(B[[1]]) + length(C[[1]][1, ])))
  type <- c(type, c(rep("ABC-Gibbs", length(A[[1]][5:50])), rep("Simple ABC", length(B[[1]])), rep("SMC-ABC", length(C[[1]][1, ]))))
}


gib_curves30 <- data.frame(
  Group = nbr,
  mu1 = par1,
  mu2 = par2,
  mu3 = par3,
  hyperparameter = hyper,
  Method = type
)



g <- ggplot(data = gib_curves30, aes(x = mu1)) + theme(panel.background = element_rect(fill = "white", colour = "grey"), legend.position = "bottom")
for (i in 1:N) {
  g <- g + geom_line(data = gib_curves30[which(gib_curves30$Group == i & gib_curves30$Method == "ABC-Gibbs"), ], stat = "density", alpha = .1)
  g <- g + geom_line(data = gib_curves30[which(gib_curves30$Group == i & gib_curves30$Method == "Simple ABC"), ], stat = "density", alpha = .1)
  g <- g + geom_line(data = gib_curves30[which(gib_curves30$Group == i & gib_curves30$Method == "SMC-ABC"), ], stat = "density", alpha = .1)
}
g <- g + facet_grid(Method ~ ., scales = "free_y") + xlab(expression(mu[1]))
g <- g + geom_line(data = data.frame(mu1 = exact[[1]][, 1], Type = rep("Exact", 1000)), stat = "density", linetype = "dashed") + scale_colour_viridis_d()
g
ggsave("mu1SMC.pdf", height = 15, width = 7, unit = "cm")

f <- ggplot(data = gib_curves30, aes(x = hyperparameter)) + theme(panel.background = element_rect(fill = "white", colour = "grey"), legend.position = "bottom")
for (i in 1:N) {
  f <- f + geom_line(data = gib_curves30[which(gib_curves30$Group == i & gib_curves30$Method == "ABC-Gibbs"), ], stat = "density", alpha = 0.2)
  f <- f + geom_line(data = gib_curves30[which(gib_curves30$Group == i & gib_curves30$Method == "Simple ABC"), ], stat = "density", alpha = 0.2)
  f <- f + geom_line(data = gib_curves30[which(gib_curves30$Group == i & gib_curves30$Method == "SMC-ABC"), ], stat = "density", alpha = 0.2)
}
f <- f + facet_grid(Method ~ ., scales = "free") + xlab(expression(alpha))
f <- f + geom_line(data = data.frame(hyperparameter = exact[[2]], Type = rep("Exact", 1000)), stat = "density", linetype = "dashed") + xlim(c(0, 5)) + scale_colour_viridis_d()
f
ggsave("hyperSMC.pdf", height = 15, width = 7, unit = "cm")

# comparaion à coût non constant

A <- gibbstot(obs, rep(0, 20), 0, sigma, vari, 50, 50, 50, 10)

hyper <- numeric()
par1 <- numeric()
par2 <- numeric()
par3 <- numeric()
nbr <- numeric()
type <- numeric()
epsfin <- matrix(ncol = 3, nrow = 10)
exact <- gibbsexact(obs, 1000, sigma, vari, 10)


N <- 10
tps <- numeric(10)
for (i in 1:N) {
  A <- gibbstot(obs, rep(0, 20), 0, sigma, vari, 30, 30, 1000, 10)
  B <- ABCsimple(obs, vari, sigma, 1000, 30000, 10)
  C <- SMCmaison(1000, 5, 30, A[[3]][length(A[[3]])], .9, obs, monprior, model)
  tps[i] <- sum(C[[7]])
  hyper <- c(hyper, A[[1]][5:50], B[[1]], C[[1]][1, ])
  par1 <- c(par1, A[[2]][1, 5:50], B[[2]][1, ], C[[1]][2, ])
  par2 <- c(par2, A[[2]][2, 5:50], B[[2]][2, ], C[[1]][3, ])
  par3 <- c(par3, A[[2]][3, 5:50], B[[2]][3, ], C[[1]][4, ])
  epsfin[1, ] <- c(A[[3]][length(A[[3]])], B[[3]], C[[3]][length(C[[3]])])
  nbr <- c(nbr, rep(i, length(A[[1]][5:50]) + length(B[[1]]) + length(C[[1]][1, ])))
  type <- c(type, c(rep("ABC-Gibbs", length(A[[1]][5:50])), rep("Simple ABC", length(B[[1]])), rep("SMC-ABC", length(C[[1]][1, ]))))
}


gib_curves30 <- data.frame(
  Group = nbr,
  mu1 = par1,
  mu2 = par2,
  mu3 = par3,
  hyperparameter = hyper,
  Method = type
)



g <- ggplot(data = gib_curves30, aes(x = mu1)) + theme(panel.background = element_rect(fill = "white", colour = "grey"), legend.position = "bottom")
for (i in 1:N) {
  g <- g + geom_line(data = gib_curves30[which(gib_curves30$Group == i & gib_curves30$Method == "ABC-Gibbs"), ], stat = "density", alpha = .2)
  g <- g + geom_line(data = gib_curves30[which(gib_curves30$Group == i & gib_curves30$Method == "Simple ABC"), ], stat = "density", alpha = .2)
  g <- g + geom_line(data = gib_curves30[which(gib_curves30$Group == i & gib_curves30$Method == "SMC-ABC"), ], stat = "density", alpha = .2)
}
g <- g + facet_grid(Method ~ ., scales = "free_y") + xlab(expression(mu[1]))
g <- g + geom_line(data = data.frame(mu1 = exact[[1]][, 1], Type = rep("Exact", 1000)), stat = "density", linetype = "dashed") + scale_colour_viridis_d()
g
ggsave("mu1SMCmieux.pdf", height = 15, width = 7, unit = "cm")

f <- ggplot(data = gib_curves30, aes(x = hyperparameter)) + theme(panel.background = element_rect(fill = "white", colour = "grey"), legend.position = "bottom")
for (i in 1:N) {
  f <- f + geom_line(data = gib_curves30[which(gib_curves30$Group == i & gib_curves30$Method == "ABC-Gibbs"), ], stat = "density", alpha = 0.2)
  f <- f + geom_line(data = gib_curves30[which(gib_curves30$Group == i & gib_curves30$Method == "Simple ABC"), ], stat = "density", alpha = 0.2)
  f <- f + geom_line(data = gib_curves30[which(gib_curves30$Group == i & gib_curves30$Method == "SMC-ABC"), ], stat = "density", alpha = 0.2)
}
f <- f + facet_grid(Method ~ ., scales = "free") + xlab(expression(alpha))
f <- f + geom_line(data = data.frame(hyperparameter = exact[[2]], Type = rep("Exact", 1000)), stat = "density", linetype = "dashed") + xlim(c(0, 5)) + scale_colour_viridis_d()
f
ggsave("hyperSMCmieux.pdf", height = 15, width = 7, unit = "cm")
