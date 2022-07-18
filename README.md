
# global-combinations

R code accompanying the paper ‘[Global combinations of expert
forecasts](https://arxiv.org/abs/2207.07318)’.

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
    ##            [,1]      [,2]      [,3]
    ## [1,] 0.28891328 0.2321148 0.1367394
    ## [2,] 0.05790316 0.1447987 0.2016698
    ## [3,] 0.13301709 0.1773916 0.1268085
    ## [4,] 0.16742777 0.2695034 0.2570763
    ## [5,] 0.35273870 0.1761914 0.2777059

``` r
# Fit the global optimal combination
combine(x, scheme = 'optimal', gamma = 1e6, lambda = 0)
```

    ## [[1]]
    ##           [,1]      [,2]      [,3]
    ## [1,] 0.2424916 0.2424911 0.2424906
    ## [2,] 0.1312026 0.1312029 0.1312034
    ## [3,] 0.1374763 0.1374767 0.1374763
    ## [4,] 0.2192649 0.2192658 0.2192654
    ## [5,] 0.2695646 0.2695636 0.2695642

``` r
# Fit the global optimal combination with grouped tasks
combine(x, scheme = 'optimal', group = c(1, 1, 2), gamma = 1e6, lambda = 0)
```

    ## [[1]]
    ##            [,1]       [,2]      [,3]
    ## [1,] 0.26967240 0.26967175 0.1367394
    ## [2,] 0.09016344 0.09016382 0.2016698
    ## [3,] 0.14355691 0.14355719 0.1268085
    ## [4,] 0.22629704 0.22629789 0.2570763
    ## [5,] 0.27031022 0.27030934 0.2777059

``` r
# Cross-validate the tuning parameters
cv.combine(x, scheme = 'optimal')
```

    ##           [,1]      [,2]      [,3]
    ## [1,] 0.2598513 0.2000267 0.2000262
    ## [2,] 0.1086179 0.1999466 0.1999468
    ## [3,] 0.1368147 0.1999427 0.1999425
    ## [4,] 0.1970347 0.2000167 0.2000166
    ## [5,] 0.2976813 0.2000672 0.2000679
