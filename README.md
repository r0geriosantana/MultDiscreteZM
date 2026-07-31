# MultDiscreteZM <img src="man/figures/logo_3d.png" align="right" height="139" />

The **MultDiscreteZM** package provides tools for simulating and fitting the Trivariate Zero-Modified COM-Poisson (ZM-CMP) model. It is designed to handle multivariate count data with zero-inflation or zero-deflation, while correctly modeling overdispersion and underdispersion through the CMP distribution.

## 📦 Installation

Currently, the **MultDiscreteZM** package is not available on CRAN, but you can install the development version directly from GitHub. 

If you do not have the `devtools` package installed on your computer, you will need to install it first:

```r
# Install devtools if necessary
# install.packages("devtools")

# Install the MultDiscreteZM package
devtools::install_github("r0geriosantana/MultDiscreteZM")
library(MultDiscreteZM) 

# 1. Define true parameters (Zero-inflation p = 0.85)
mu1 <- 3.7; mu2 <- 5.8; mu3 <- 6.8
kappa1 <- -0.25; kappa2 <- -0.30; kappa3 <- -0.21
p <- 0.85

# 2. Simulate a single dataset (n = 100)
set.seed(123)
Y <- r3zmcmpIu(n = 100, mu1, mu2, mu3, kappa1, kappa2, kappa3, p)

# 3. Fit the Trivariate ZMCMP Model
fit <- TrizmcmpIu(Y1 = Y[,1], Y2 = Y[,2], Y3 = Y[,3])

# The function will output the estimated parameters and standard errors!
## 📊 Example 2: Regression Model with Covariates

# 1. Define sample size and simulate a continuous covariate 'x'
set.seed(131)
n <- 100
x <- round(runif(n, 0.7, 0.9), 2)

# 2. Define parameters for the regression (beta, dispersion kappa, and zero-modification w0)
beta1 <- c(0.83, -0.79); kappa1 <- -1.09
beta2 <- c(0.80, -0.76); kappa2 <- -1.12
beta3 <- c(0.74, -0.70); kappa3 <- -1.13
w0 <- 0.95

# 3. Simulate multivariate regression data
Y <- r3zmcmpuI_reg0(x, beta1, beta2, beta3, kappa1, kappa2, kappa3, w0)

# 4. Fit the Trivariate ZMCMP Regression Model
fit_reg <- reg3zmcmpuI(Y1 = Y[,1,], Y2 = Y[,2,], Y3 = Y[,3,], x = x)
