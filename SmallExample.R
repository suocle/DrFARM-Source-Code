n <- 500
p <- 10
q <- 5

i <- 1
set.seed(i)

## generate X matrix
X <- matrix(rnorm(n*p), n, p)

## generate coefficient
Theta.t <- matrix(0, p, q)

idx <- matrix(c(rep(3, 3), rep(8, 2), rep(10, 4),
                1, 2, 4, 3, 5, 1, 2, 4, 5), 9, 2)

for (i in 1:9) {
  rand <- runif(1, min = 0, max = 1) > 0.5
  rand[rand == 0] <- -1
  coef <- runif(1, min = 1, max = 6) * rand
  Theta.t[idx[i,1], idx[i,2]] <- coef
}

k <- 2
tau = 11.49
B.star <- matrix(runif(q*k, min = 0, max = tau), q, k) #tau = 11.49 for this simple example
BBT <- B.star %*% t(B.star)
eigen.res <- eigen(BBT)
U <- eigen.res$vectors
U <- U[,1:k]
v <- eigen.res$values
v <- v[1:k]
sqrt.V <- diag(sqrt(v))
B <- U %*% sqrt.V

Z <- matrix(rnorm(n*k), n, k)
r <- mean(diag(cov(X %*% Theta.t))) / mean(diag(cov(Z %*% t(B))))

#Signal to signal to noise ratio (SSNR) set to 1:3:5
sigma <- sqrt(5 * mean(diag(cov(X %*% Theta.t))))
E <- matrix(rnorm(n*q, sd = sigma), n, q)
Y <- X %*% Theta.t + Z %*% t(B) + E

rm(list=setdiff(ls(), c("X", "Y", "Theta.t")))
