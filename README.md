
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
set.seed(2022)
n <- 30 # Number of forecast errors
p <- 5 # Number of forecasters
m <- 3 # Number of tasks
x <- replicate(m, matrix(rnorm(n * p), n, p), simplify = F)

# Fit the local optimal combination
combine(x, scheme = 'optimal', gamma = 0, lambda = 0)
```

    ## [[1]]
    ##           [,1]      [,2]       [,3]
    ## [1,] 0.1333955 0.1578834 0.17471355
    ## [2,] 0.3427564 0.2478892 0.26697511
    ## [3,] 0.1569549 0.1896682 0.38160184
    ## [4,] 0.2631108 0.1629942 0.05353802
    ## [5,] 0.1037825 0.2415651 0.12317148

``` r
# Fit the global optimal combination
combine(x, scheme = 'optimal', gamma = 1e6, lambda = 0)
```

    ## [[1]]
    ##           [,1]      [,2]      [,3]
    ## [1,] 0.1573271 0.1573271 0.1573272
    ## [2,] 0.2694190 0.2694186 0.2694188
    ## [3,] 0.2307956 0.2307959 0.2307967
    ## [4,] 0.1737718 0.1737713 0.1737708
    ## [5,] 0.1686866 0.1686871 0.1686865

``` r
# Fit the global optimal combination with grouped tasks
combine(x, scheme = 'optimal', group = c(1, 1, 2), gamma = 1e6, lambda = 0)
```

    ## [[1]]
    ##           [,1]      [,2]       [,3]
    ## [1,] 0.1462684 0.1462685 0.17471355
    ## [2,] 0.2744725 0.2744721 0.26697511
    ## [3,] 0.1816342 0.1816343 0.38160184
    ## [4,] 0.2216673 0.2216668 0.05353802
    ## [5,] 0.1759576 0.1759582 0.12317148

``` r
# Cross-validate the tuning parameters
cv.combine(x, scheme = 'optimal')
```

    ##           [,1]      [,2]      [,3]
    ## [1,] 0.1999658 0.1999658 0.1909410
    ## [2,] 0.2000542 0.2000540 0.2373517
    ## [3,] 0.2000331 0.2000333 0.3093194
    ## [4,] 0.1999739 0.1999735 0.1114921
    ## [5,] 0.1999730 0.1999734 0.1508959
