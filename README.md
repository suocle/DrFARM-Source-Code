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

which contains a small toy example with sample size `n = 500`, `p = 10` variants and `q = 5` traits, with 3 pleiotropic variants (variant #3, #8 and #10). For simplicity, most objects from the R environment are removed, except the main input required: `X` (`n x p` variants matrix) and `Y` (`n x q` trait matrix), as well as the ground truth `p x q` coefficient matrix `Theta.t`.
