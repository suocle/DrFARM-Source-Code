# DrFARM-Source-Code

### Overview
This repository provides a demonstration on how to use the DrFARM `R` source code.

# System Requirements

## Software Requirements

### OS Requirements

The source code has been tested on Microsoft's Windows 10 operating system and Linux (Ubuntu 18.04). The source code should be compatible with Windows, Mac, and Linux operating systems.

Before using the DrFARM source code, users should have `R` version 4.0.4 or higher, and several packages installed.

### Installation  

First, we need to install `glmnet`, `glasso` and `psych`:  

    install.packages(c("glmnet", "glasso", "psych"))
    
which should install within a couple of minutes on a standard machine.
   
We also need to install a modified version of `remMap` (which is no longer on CRAN) available on this repository:

    install.packages("remMap_0.2-0.tar.gz", repos = NULL)
    
If you are interested in replicating simulation II, you will also need to install `mvnfast`:  

    install.packages("mvnfast")

# Demo

We first load all the source code dependencies:

```
library(glmnet)
library(glasso)
library(psych)
library(remMap)
#library(mvnfast) #for replicating simulation II only
```

and the source code containing all the main functions:

```
source("DrFARM.R")
```

Then, we source 

```
source("SmallExample.R")
```

which contains a small toy example with sample size `n = 500`, `p = 10` variants and `q = 5` traits, with 3 pleiotropic variants (variant #3, #8 and #10). For simplicity, most objects from the R environment are removed, except the main input required: `X` (`n x p` variants matrix) and `Y` (`n x q` trait matrix), as well as the ground truth `p x q` coefficient matrix `Theta.t`. Notice that there are `k = 2` underlying latent factors that contribute to the 5 traits.

We estimate the initial value `Theta0` using remMap<sup>[1]</sup>:
```
remMap.res <- remMap.whole(X, Y)
Theta0 <- remMap.res$Theta0
```

Then, we estimate the precision matrix using `precM`:
```
precM <- precM(X)
```
By default, `precM` estimates the precision matrix using glasso (the recommended approach in our paper).

In this example, we assume the number of latent factors is known (`k = 2`). Using `DrFARM.whole`:
```
k <- 2
DrFARM.res <- DrFARM.whole(X, Y, Theta0, precM, k, 
                           remMap.res$lambda1.opt, 
                           remMap.res$lambda2.opt)
Theta <- DrFARM.res$Theta;
B <- DrFARM.res$B; 
E.Z <- DrFARM.res$E.Z;
```
we obtain the estimated `q x p` sparse coefficient matrix `Theta`, `q x k` loading matrix `B` and `n x k` expected latent factors. These output are essential for the final step of calculating the entrywise *p*-values as well as the pleiotropic *p*-values.

The `q x p` entrywise (`pval1`) and length `p` pleiotropic (`pval2`) *p*-values can simply be obtained using:
```
pval1 <- entry.pvalue(X, Y, Theta, B, E.Z, precM)
pval2 <- pleio.pvalue(X, Y, Theta, B, E.Z, precM)
```

# Output
```
pval1
          [,1]       [,2]         [,3]      [,4]        [,5]      [,6]      [,7]         [,8]      [,9]        [,10]
[1,] 0.7490623 0.88137417 5.824796e-19 0.3932252 0.070366934 0.3090272 0.3464442 6.690187e-01 0.6855273 5.990873e-13
[2,] 0.0783111 0.54949451 2.550019e-07 0.7744600 0.645618001 0.1126868 0.7932174 5.272138e-01 0.5744373 6.029076e-26
[3,] 0.8307014 0.14338934 1.088501e-01 0.4906801 0.898785814 0.3863579 0.2613619 7.648873e-25 0.5073458 3.025948e-01
[4,] 0.7089323 0.09199648 1.767313e-07 0.9791283 0.284369858 0.8665795 0.2215291 1.049544e-01 0.6258750 1.256087e-34
[5,] 0.4909801 0.07027779 6.875141e-01 0.3353919 0.007001236 0.1024763 0.9474172 1.175925e-02 0.6030415 5.418847e-51

pval2
[1] 9.306550e-01 3.875994e-01 1.949086e-16 1.992648e-01 6.721478e-02 5.323137e-01 5.221915e-01 1.949086e-16 7.943848e-01 6.496883e-17
```

# Simulation I replication
In simulation I of our manuscript<sup>[2]</sup> (which assumes the individuals to be unrelated), interested users can easily re-generate a replicate of the debiasing-based results by changing several parameters from the above toy example code, namely, `n`, `p`, `q`, `k`, `tau`, replace line 14-22 of `SmallExample.R`:

```
idx <- matrix(c(rep(3, 3), rep(8, 2), rep(10, 4),
                1, 2, 4, 3, 5, 1, 2, 4, 5), 9, 2)

for (i in 1:9) {
  rand <- runif(1, min = 0, max = 1) > 0.5
  rand[rand == 0] <- -1
  coef <- runif(1, min = 1, max = 6) * rand
  Theta.t[idx[i,1], idx[i,2]] <- coef
}
```

by

```
#Number of pleiotropic variants
pv <- ##

#Number of signals
signal <- ##

#Index of pleiotropic variants
pv.idx <- sort(sample(1:p, pv))

#Number of traits regulated per pleiotropic variant
t.per.pv <- rmultinom(1, signal, rep(1/pv, pv))
for (i in 1:pv) {
  trait.idx <- sort(sample(1:q, t.per.pv[i]))
  rand <- runif(t.per.pv[i], min = 0, max = 1) > 0.5
  rand[rand == 0] <- -1
  coef <- runif(t.per.pv[i], min = 1, max = 1.5) * rand
  Theta.t[pv.idx[i], trait.idx] <- coef
}
```

and varying the choice of precision matrix in `precM`: Glasso<sup>[3]</sup> (`method = "glasso"`), nodewise lasso<sup>[4]</sup> (`method = "NL"`) and quadratic optimization<sup>[5]</sup> (`method = "QO"`). For scenario I, use `n = 1000`, `p = 2000`, `q = 500`, `k = 5`, `pv = 300`, `signal = 3000` and `tau = 4.151`; for scenario II, use `n = 2000`, `p = 5000`, `q = 1000`, `k = 10`, `pv = 750`, `signal = 7500` and `tau = 3.276`. However, given the computational complexity of the problem, it is recommended to use the parallizable counterparts (`remMap.one` and `DrFARM.one`) of the key functions (`remMap.whole` and `DrFARM.whole`), and use a computer cluster (e.g., array jobs) in order to achieve a faster result (see below).

# Helper functions used in parallization
Unlike `remMap.whole` and `DrFARM.whole`, which are all-in-functions, we need to manually generate the tuning parameter grid. This can be done via the use of `remMap.grid` and `DrFARM.grid`. By default, `remMap.grid` generates a 10 x 10 grid (in a similar manner to how `glmnet` generate a sequence of 100 tuning parameter values) and `DrFARM.grid` generates a 5 x 5 grid (described in our manuscripts appendix<sup>[2]</sup>). Continuning our toy example from above, below shows an example of remMap using `i = 36` (36th row of `remMap.lambda.grid`, which also happens to the optimal pairs of tuning parameters).

```
remMap.lambda.grid <- remMap.grid(X, Y)

i <- 36

#Initial value candidate
Theta0.cand <- remMap.one(X, Y, lambda1 = remMap.lambda.grid[i,1], lambda2 = remMap.lambda.grid[i,2])
```

Assuming we obtained all the 100 (10 x 10) initial value candidates using `remMap.one`, one possible way to select the remMap initial value is choosing the one that minimizes the extended Bayesian information criterion<sup>[6]</sup> (EBIC). The calculation of EBIC for remMap can be done using

```
EBIC <- remMap.EBIC(X, Y, Theta0.cand)
```

By default, remMap.EBIC uses `gamma = 1` (`gamma` is a hyperparameter ranging from 0 to 1 that controls how much EBIC prefers simpler models) which yields the highest possible sparsity among the choices of `gamma`. If BIC (a special case is EBIC when `gamma = 0`) is desired, this can also be calculated by specifying the option `gamma = 0`:
```
BIC <- remMap.EBIC(X, Y, Theta0.cand, gamma = 0)
```

Next, assuming both the initial value `Theta0` and precision matrix `precM` are available. We again generate a tuning parameter grid based on the tuning parameters of `Theta0`:

```
DrFARM.lambda.grid <- DrFARM.grid(X, Y, Theta0, precM, k = 2, lambda1.opt = remMap.lambda.grid[i,1], lambda2.opt = remMap.lambda.grid[i,2])
```

Notice that the tuning parameter grid generated using our method depends on the number of latent factors `k`, which is assumed to be fixed in both the toy example and simulation studies. Then, `DrFARM.one` can be similarly applied. Below shows an example using `i = 19` (19th row of `DrFARM.lambda.grid`, which also happens to the optimal pairs of tuning parameters).


```
i <- 19

#DrFARM estimates candidate
DrFARM.one.res <- DrFARM.one(X, Y, Theta0, precM, k = 2, lambda1 = DrFARM.lambda.grid[i,1], lambda2 = DrFARM.lambda.grid[i,2])
```

Finally, assume we obtained all the 25 (5 x 5) DrFARM sparse solution candidates using `DrFARM.one`, we can again select the sparse estimate that minimizes EBIC. The calculation of EBIC for DrFARM is done by
```
Theta <- DrFARM.one.res$Theta;
B <- DrFARM.one.res$B; 
E.Z <- DrFARM.one.res$E.Z;
diag.Psi <- DrFARM.one.res$diag.Psi;
EBIC <- DrFARM.EBIC(X, Y, Theta, B, E.Z, diag.Psi)
```

# Simulation II replication
In simulation II of our manuscript<sup>[2]</sup>, the individuals are assumed to have third-degree relatedness on average. Assuming the `mvnfast` library has been loaded, in order to generate a replicate of data for simulation II, all we need is to replace line 36 of `SmallExample.R`:

```
Z <- matrix(rnorm(n*k), n, k)
```

by

```
R <- matrix(rbinom(n*p, size = 1, prob = 0.25), n, p, byrow = TRUE)
K <- cov2cor(R %*% t(R))
Z <- t(rmvn(k, rep(0, n), K))
```

and line 44:

```
rm(list=setdiff(ls(), c("X", "Y", "Theta.t")))
```

by

```
rm(list=setdiff(ls(), c("X", "Y", "Theta.t", "K)))
```

and changing several parameters just as in simulation I. For scenario I, use `n = 1000`, `p = 2000`, `q = 500`, `k = 5`, `pv = 300`, `signal = 3000` and `tau = 2.891`; for scenario II, use `n = 2000`, `p = 5000`, `q = 1000`, `k = 10`, `pv = 750`, `signal = 7500` and `tau = 2.281`. Running DrFARM on data with kinship only requires specifying one additional argument `K = K` (by default `K = NULL`) in `DrFARM.grid`, `DrFARM.one`, `DrFARM.EBIC` and `DrFARM.all`. For example,

```
DrFARM.res <- DrFARM.whole(X, Y, Theta0, precM, k, 
                           remMap.res$lambda1.opt, 
                           remMap.res$lambda2.opt,
                           K = K)
```

### References

[1] Peng, J., Zhu, J., Bergamaschi, A., Han, W., Noh, D. Y., Pollack, J. R., & Wang, P. (2010). Regularized multivariate regression for identifying master predictors with application to integrative genomics study of breast cancer. The annals of applied statistics, 4(1), 53.

[2] Chan, L. S., Li, G., Fauman, E. B., Laakso, M., Boehnke, M., & Song, P. X. (2022). DrFARM: Identification and inference for pleiotropic gene in GWAS. bioRxiv.

[3] Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse covariance estimation with the graphical lasso. Biostatistics, 9(3), 432-441.

[4] Van de Geer, S., B??hlmann, P., Ritov, Y. A., & Dezeure, R. (2014). On asymptotically optimal confidence regions and tests for high-dimensional models. The Annals of Statistics, 42(3), 1166-1202.

[5] Javanmard, A., & Montanari, A. (2014). Confidence intervals and hypothesis testing for high-dimensional regression. The Journal of Machine Learning Research, 15(1), 2869-2909.

[6] Chen, J., & Chen, Z. (2008). Extended Bayesian information criteria for model selection with large model spaces. Biometrika, 95(3), 759-771.
