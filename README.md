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

# Demo

We first load all the source code dependencies:

```
library(glmnet)
library(glasso)
library(psych)
library(remMap)
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

We estimate the initial value `Theta0` using remMap:
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
In simulation I of our manuscript (which assumes the individuals to be unrelated), the debiasing-based results can be replicated by changing several parameters, namely, `n`, `p`, `q`, `k`, `tau` and varying the choice of precision matrix in `precM`: glasso (`method = "glasso"`), nodewise lasso (`method = "NL"`) and quadratic optimization (`method = "QO"`).
