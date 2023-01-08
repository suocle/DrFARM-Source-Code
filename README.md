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

## Functions

For interactive demos of the functions, please check out the vignettes built into the package. They can be accessed as follows:

```
require(lolR)
vignette('lol')
vignette('pca')
vignette('cpca')
vignette('lrcca')
vignette('mdp')
vignette('xval')
vignette('qoq')
vignette('simulations')
vignette('nearestCentroid')
```

## Extending the lolR Package

The lolR package makes many useful resources available (such as embedding and cross-validation) for simple extension. 

To extend the lolR package, check out the vignettes:

```
require(lolR)
vignette('extend_embedding')
vignette('extend_classification')
```

# Results

In this [benchmark comparison](http://docs.neurodata.io/lol/lol-paper/figures/real_data.html), we show that LOL does better than all linear embedding techniques in supervised HDLSS settings when dimensionality is high (d > 100, ntrain <= d) on 20 benchmark problems from the [UCI](https://archive.ics.uci.edu/ml/index.php) and [PMLB](https://github.com/EpistasisLab/penn-ml-benchmarks) datasets. LOL provides a good tradeoff between maintaining the class conditional difference (good misclassification rate) in a small number of dimensions (low number of embedding dimensions).
