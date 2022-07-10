
# global-combinations

R code accompanying the paper ‘Global combinations of expert forecasts’.

You should install Gurobi and the associated R package `gurobi` before
using this code. Gurobi is available for free under academic license at
<https://www.gurobi.com/>.

## Usage

The `combine()` function fits forecast combination weights for a
sequence of tuning parameters. The tuning parameter `gamma` controls the
level of globalisation and `lambda` controls the level of shrinkage. The
`cv.combine` function chooses good values of the parameters
automatically using leave-one-out cross-validation.

``` r
source('combination-functions.R')

# Generate forecast errors for multiple forecasters across several tasks
n <- 30 # Number of forecast errors
p <- 5 # Number of forecasters
m <- 3 # Number of tasks
x <- replicate(m, matrix(rnorm(n * p), n, p), simplify = F)

# Fit the local optimal combination
combine(x, scheme = 'optimal', gamma = 0, lambda = 0)
```

    ## [[1]]
    ##            [,1]      [,2]       [,3]
    ## [1,] 0.13412812 0.2289420 0.33878341
    ## [2,] 0.02864837 0.2520558 0.29258236
    ## [3,] 0.35281229 0.2646848 0.10041654
    ## [4,] 0.26603953 0.1110819 0.24978408
    ## [5,] 0.21837168 0.1432355 0.01843361

``` r
# Fit the global optimal combination
combine(x, scheme = 'optimal', gamma = 1e6, lambda = 0)
```

    ## [[1]]
    ##           [,1]      [,2]      [,3]
    ## [1,] 0.2165707 0.2165714 0.2165718
    ## [2,] 0.2021712 0.2021719 0.2021721
    ## [3,] 0.2712285 0.2712275 0.2712270
    ## [4,] 0.1966654 0.1966650 0.1966655
    ## [5,] 0.1133643 0.1133643 0.1133636

``` r
# Fit the global optimal combination with grouped tasks
combine(x, scheme = 'optimal', group = c(1, 1, 2), gamma = 1e6, lambda = 0)
```

    ## [[1]]
    ##           [,1]      [,2]       [,3]
    ## [1,] 0.1850168 0.1850174 0.33878341
    ## [2,] 0.1527581 0.1527589 0.29258236
    ## [3,] 0.3175761 0.3175754 0.10041654
    ## [4,] 0.1898926 0.1898919 0.24978408
    ## [5,] 0.1547565 0.1547565 0.01843361

``` r
# Cross-validate the tuning parameters
cv.combine(x, scheme = 'optimal')
```

    ##           [,1]      [,2]      [,3]
    ## [1,] 0.1586639 0.2000175 0.2345242
    ## [2,] 0.1049071 0.1999988 0.2228004
    ## [3,] 0.3128793 0.2000590 0.1820349
    ## [4,] 0.2321819 0.2000052 0.2113089
    ## [5,] 0.1913679 0.1999195 0.1493315
