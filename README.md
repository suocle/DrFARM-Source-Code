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
