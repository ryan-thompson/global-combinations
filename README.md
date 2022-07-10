
# global-combinations

Code accompanying the paper ‘Global combinations of expert forecasts’.

You should install Gurobi and the associated R package gurobi before
using this code. Gurobi is available for free under academic license at
<https://www.gurobi.com/>.

## Usage

The `combine()` function fits forecast combination weights for a
sequence of tuning parameters. The tuning parameter `gamma` controls the
level of globalisation and `lambda` controls the level of shrinkage. The
`cv.combine` functions chooses the best values of the parameters
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
    ##           [,1]      [,2]      [,3]
    ## [1,] 0.0909584 0.2702194 0.1364487
    ## [2,] 0.2738042 0.1494573 0.1877587
    ## [3,] 0.2652071 0.1053267 0.2954725
    ## [4,] 0.1998196 0.2049365 0.1845877
    ## [5,] 0.1702108 0.2700600 0.1957324

``` r
# Fit the global optimal combination
combine(x, scheme = 'optimal', gamma = 1e6, lambda = 0)
```

    ## [[1]]
    ##           [,1]      [,2]      [,3]
    ## [1,] 0.1640282 0.1640289 0.1640282
    ## [2,] 0.1900090 0.1900085 0.1900087
    ## [3,] 0.2355067 0.2355058 0.2355070
    ## [4,] 0.1938856 0.1938859 0.1938858
    ## [5,] 0.2165705 0.2165708 0.2165703

``` r
# Fit the global optimal combination with grouped tasks
combine(x, scheme = 'optimal', gamma = 1e6, lambda = 0, group = c(1, 1, 2))
```

    ## [[1]]
    ##           [,1]      [,2]      [,3]
    ## [1,] 0.1925696 0.1925703 0.1364487
    ## [2,] 0.2024248 0.2024244 0.1877587
    ## [3,] 0.1860809 0.1860800 0.2954725
    ## [4,] 0.1805066 0.1805069 0.1845877
    ## [5,] 0.2384181 0.2384184 0.1957324

``` r
# Cross-validate the globalisation parameter
cv.combine(x, scheme = 'optimal', lambda = 0)
```

    ##           [,1]      [,2]      [,3]
    ## [1,] 0.1638055 0.1988963 0.1638452
    ## [2,] 0.1903571 0.1798172 0.1899865
    ## [3,] 0.2355939 0.1849486 0.2359076
    ## [4,] 0.1936921 0.1994283 0.1938871
    ## [5,] 0.2165514 0.2369096 0.2163735
