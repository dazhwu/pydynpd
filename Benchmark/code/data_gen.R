#code is adapted from 
#Phillips PC, Han C (2019). “Chapter 5 - Dynamic panel GMM using R.” In HD Vinod, C Rao
#(eds.), Conceptual Econometrics Using R, volume 41 of Handbook of Statistics, pp. 119 –
#144. Elsevier. doi:https://doi.org/10.1016/bs.host.2019.01.002.

library(Matrix)
library(data.table)


nsize = 1000
tsize = 10
burn = 20
t.all = tsize + burn

beta1 = .25
beta2 = .1
gamma0 = 5
gamma1 = -.1
gamma2 = .1
gamma3 = 3
sigma = list(alpha = 5, u = 1, x = 1, z = 1)
mu0 = 1
mu1 = .005
rho = list(z = 1, x = 1)

set.seed(1)
x = sigma$x * matrix(rnorm(nsize * t.all), nsize, t.all)
for (j in 2:ncol(x)) x[, j] = rho$x * x[, j - 1] + x[, j]
if (rho$x == 1) x = x - x[, burn - 1]

trend = mu0 + mu1 * seq(-burn + 1, tsize)
z0 = stats::filter(rnorm(t.all), rho$z, method = "recursive")
if (rho$z == 1) z0 = z0 - z0[burn - 1]
z = trend + z0

alpha = sigma$alpha * rnorm(nsize)
u = sigma$u * matrix(rnorm(nsize * t.all), nsize, t.all)

# With these components in hand, we recursively generate yit as follows:
xbar = colMeans(x)
y = matrix(NA, nsize, t.all)
y[, 1] = alpha + u[, 1]
for (j in 2:ncol(y)) {
  lambda = gamma0 +     gamma1 * mean(y[, j - 1]) +     gamma2 * xbar[j - 1] +     gamma3 * z[j - 1]
  y[, j] = alpha +     beta1 * y[, j - 1] +     beta2 * x[, j - 1] +     lambda +     u[, j]
}

y = y[, burn:t.all]  # y is the Nx(T+1) matrix    t=0...T
x = x[, burn:t.all]
z = z[burn:t.all]

w=data.frame(id=as.vector(row(y)), year=as.vector(col(y))-1, y=as.vector(y), x=as.vector(x))
write.csv(w, "test_data.csv")