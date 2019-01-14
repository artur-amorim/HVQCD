[![license](https://img.shields.io/github/license/artur-amorim/HVQCD.svg)]()
[![GitHub language count](https://img.shields.io/github/languages/count/artur-amorim/HVQCD.svg)]()
[![GitHub issues](https://img.shields.io/github/issues/artur-amorim/HVQCD.svg)]()
# HVQCD
An R package that computes the numerical solution of the holographic models of QCD in the Veneziano limit by Jarvinen and Kiritsis. R is used as an interface between the C++ code and the user.

# Dependencies

## Boost
This packages requires the [Boost](https://www.boost.org/) libraries on your computer. They were used to compute the solutions of the equations of motion of the background fields. Look in the official website for more details on how to install it.

## Rcpp
This library, among other things, allows you to call C++ functions in R.

# Install
Type the following in R in order to install this tool
```r
# install devtools if you do not have it
install.packages('devtools')
# install the tool
devtools::install_github('artur-amorim/HVQCD')
```

### Usage
Here's how to use it in R

```r
library(HVQCD)
x <- 1.0
t0 <- 1.0
# computs the profile of the warp factor A, of the dilaton and of the tachyon
s <- solveHVQCD(x, t0)
# s$z: The z values
# s$A: The corresponding A values
# s$lambda : The corresponding lambda = exp(dilaton) values
# s$tau : The corresponding values of the tachyon
# s$mq : Value of the quark mass associated with the parameters x and 
# plot of A as a function of z
plot(s$z, sol$A)
```
# References
- Holographic Models for QCD in the Veneziano Limit [arXiv:1112.1261v2](https://arxiv.org/abs/1112.1261)