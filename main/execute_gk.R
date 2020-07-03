
library(gk)
library(smfsb)
library(ggplot2)
library(mvtnorm)

source("functions_gk.R")

# simulation de données


hyper <- runif(1, -10, 10)
As <- rnorm(50, hyper, 1)
B <- runif(1, 0, 1)
g <- runif(1, 0, 1)
k <- runif(1, 0, 1)

X <- matrix(NA, ncol = 50, nrow = 20)
for (i in 1:50) {
  X[, i] <- rgk(20, As[i], B, g, k)
}


Statstar <- numeric()
for (i in 1:50) {
  Statstar <- c(Statstar, quantile(X[, i], (0:8) / 8))
}

# définition du modèle pour SMCmaison
rdist <- function(par, Statstar) {
  stat <- numeric()
  for (i in 1:50) {
    stat <- c(stat, quantile(rgk(20, par[1 + i], par[52], par[53], par[54]), (0:8) / 8))
  }
  return(sum(abs(stat - Statstar)))
}

rprior <- function() {
  P <- numeric(54)
  P[1] <- runif(1, -10, 10)
  P[52:54] <- runif(3, 0, 1)
  P[2:51] <- rnorm(50, P[1], 1)
  return(P)
}

dprior <- function(P, ...) {
  dunif(P[1], -10, 10, ...) + sum(dunif(P[52:54], 0, 1, ...)) + sum(dnorm(P[2:51], P[1], 1, ...))
}


model <- rdist

monprior <- list(density = dprior, simu = rprior)




# comparaison

outSMC <- SMCmaison(1000, 5, 500, 0, .9, Statstar, monprior, model)
outVan <- abcvan(1000, 200 * 1000, Statstar)
outgib <- gibbs(c(100, 50, 50, 50, 50), X, 1000)


# plots


Dat <- data.frame(
  value = c(
    outSMC[[1]][2, ], outVan[, 2], outgib[, 2],
    outSMC[[1]][3, ], outVan[, 3], outgib[, 3],
    outSMC[[1]][4, ], outVan[, 4], outgib[, 4],
    outSMC[[1]][5, ], outVan[, 5], outgib[, 5]
  ),
  Method = rep(c(
    rep("ABC-SMC", 1000), rep("vanilla ABC", 1000),
    rep("ABC Gibbs", 1000)
  ), 4),
  Parameter = c(rep("mu1", 3000), rep("mu2", 3000), rep("mu3", 3000), rep("mu4", 3000))
)

Dat$Parameter <- factor(Dat$Parameter,
  levels = c("mu1", "mu2", "mu3", "mu4"),
  labels = c(
    "mu1" = expression(mu[1]),
    "mu2" = expression(mu[2]),
    "mu3" = expression(mu[3]),
    "mu4" = expression(mu[4])
  )
)

theta_star <- data.frame(
  value = As[1:4],
  Parameter = levels(Dat$Parameter)
)

f <- ggplot(data = Dat) + stat_density(aes(x = value, linetype = Method),
  geom = "line", position = "identity"
) +
  geom_vline(data = theta_star, aes(xintercept = value), colour = "red") +
  theme(panel.background = element_rect(fill = "white", colour = "grey"), legend.position = "bottom") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
f <- f + coord_cartesian(xlim = c(0, 4), expand = TRUE, default = FALSE, clip = "on")
# geom_density(data=Dat,aes(x=value,color=Method))

f <- f + facet_grid(. ~ Parameter, labeller = label_parsed) + scale_colour_viridis_d()
f
ggsave("4parbis_nouv.pdf", height = 7, width = 15, units = "cm")

# plots pour les hyperparamètres

Dat <- data.frame(
  value = c(
    outSMC[[1]][1, ], outVan[, 1], outgib[, 1],
    outSMC[[1]][52, ], outVan[, 52], outgib[, 52],
    outSMC[[1]][53, ], outVan[, 53], outgib[, 53],
    outSMC[[1]][54, ], outVan[, 54], outgib[, 54]
  ),
  Method = rep(c(
    rep("ABC-SMC", 1000), rep("vanilla ABC", 1000),
    rep("ABC Gibbs", 1000)
  ), 4),
  Parameter = c(rep("hyperparameter", 3000), rep("B", 3000), rep("g", 3000), rep("k", 3000))
)

Dat$Parameter <- factor(Dat$Parameter,
  levels = c("hyperparameter", "B", "g", "k"),
  labels = c(
    "hyperparameter" = expression(alpha),
    "B" = expression(B),
    "g" = expression(g),
    "k" = expression(k)
  )
)

theta_star <- data.frame(
  value = c(hyper, B, g, k),
  Parameter = levels(Dat$Parameter)
)

f <- ggplot(data = Dat) + stat_density(aes(x = value, linetype = Method),
  geom = "line", position = "identity"
) +
  geom_vline(data = theta_star, aes(xintercept = value), colour = "red") +
  theme(panel.background = element_rect(fill = "white", colour = "grey"), legend.position = "bottom") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
# geom_density(data=Dat,aes(x=value,color=Method))
f <- f + coord_cartesian(ylim = c(0, 3), expand = TRUE, default = FALSE, clip = "on")
f <- f + facet_grid(. ~ Parameter, labeller = label_parsed, scale = "free")
f

ggsave("autresparbis_nouv.pdf", height = 7, width = 15, units = "cm")
