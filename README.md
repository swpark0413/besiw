# besiw

### Overview
-------


The `bisew` package allows you to estimate eigenvalues and eigenvectors and fit a principal component regression (PCR) using a Bayesian approach  under the assumption of a spiked high-dimensional covariance matrix. A key advantage of the proposed models is that they can provide uncertainty estimates. The algorithm is also straightforward to implement, as it requires only inverse-Wishart sampling and spectral decomposition.


### Installation
-------

You can install the package in two ways: 

#### Option 1: Install the most up-to-date version via `devtools`.

To install the most up-to-date version from GitHub, use the following commands:


``` R
if (!("devtools" %in% installed.packages()[,"Package"])) {
    install.packages("devtools")
}
devtools::install_github("swpark0413/besiw")
```

#### Option 2: Install from a bundled package.

First, download a bundled package from the [Releases](https://github.com/swpark0413/besiw/releases/) page. Then, install it using the command below:

``` R
# Replace "~/path/bisew_0.0.0.tar.gz" with the path to your downloaded file
install.packages("~/path/bisew_0.0.0.tar.gz", type = "source", repos = NULL)
```


### Example
-------

Here is an example of how to use the `bisew` package:

``` r
library(bisew)

# generate a spiked covariance matrix:
n <- 50
p <- 100
K <- 3

leading <- c(50,20,10)
remaining <- rep(1, p - K)
Sigma0 <- diag(c(leading, remaining), p)

# generate data
set.seed(413)
X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma0)
 
# Compute both posterior means and credible intervals for eigenvalues and eigenvectors of a covariance matrix:
H <- 4 * diag(rep(1, p))
iter <- 100000
burnin <- floor(iter / 2)
thin <- 100
res <- besiw::eigenGSIW(X = X, k = K, H = H, nmcmc = iter, nburn = burnin, nthin = thin)
est <- besiw::estimate(res)
```



### Citation
-------

If you use the bisec package in your research, please cite the following paper:

- Seongmin Kim, Kwangmin Lee, Sewon Park, and Jaeyong Lee.
  Eigenstructure inference for high-dimensional covariance with generalized shrinkage inverse-Wishart prior.
  arXiv preprint arXiv:2412.10753.

<!-- BibTeX citation:
``` bibtex
@Article{ZhangRD2022gps,
  author        = {Zhang, Ruda and Mak, Simon and Dunson, David},
  title         = {Gaussian Process Subspace Prediction for Model Reduction},
  journal       = {SIAM Journal on Scientific Computing},
  year          = {2022},
  volume        = {44},
  number        = {3},
  pages         = {A1428-A1449},
  doi           = {10.1137/21M1432739},
}
``` -->


### License
-------

The package is distributed under the GPL-3.0 license. See the [`LICENSE`](LICENSE) file for more details

