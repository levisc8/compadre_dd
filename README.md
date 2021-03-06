
Attempting to create a framework for density dependence in Comapdre
-------------------------------------------------------------------

This is pretty experimental, but the idea is to store matrix elements (which may be constants *or* expressions) in a long data format, then use a couple functions to take a user- or database supplied population vector and generate a density dependent matrix. Outputs are then generated via iteration.

These functions depend on `rlang` to work. `rlang` itself has no further dependencies, so this would be a fairly lightweight addition to `Rcompadre`/`popdemo`. Unfortunately, they do make use of `env_bind_lazy()` which is currently listed as experimental in the `rlang` [lifecycle](https://rlang.r-lib.org/reference/lifecycle.html). If this is dropped in subsequent versions, we'll need to implement the delayed assignment manually.

### Density-dependent matrices

These are now implemented. `iterate_dd_mat()` can handle both user-supplied and data base matrices.

`CompadreDDM` matrices do not look like other matrices stored in Compadre. They are lists with 2 elements: a `data_list` which contains values for each parameter and a `mat_exprs` list which contains expressions to calculate density dependent vital rates (e.g. survival, growth, reproduction). Each of these can accept any number of named values and has one additional required argument. `data_list` requires a named vector called `initial_population_vector` and the `mat_exprs` list requiress an expression for the matrix (`mat_expr`). Hopefully, the example below clarifies how these functions work.

*To-do: incorporate `stringToMatrix` from `Rcompadre` so entering the matrices is easier on Compadrinos*

``` r
# This is not yet part of the Rcompadre package, so you'll need source()
# the functions

source('R/functions.R')

# create expressions for each vital rate in the matrix using make_mat_exprs()

exprs <- make_mat_exprs(
  s_2 = 1/(1 + exp(bs2_2 * u_i + bs2_1 * t_i + bs2_0)),
  s_3 = exp(bs3_1 * log(r + 1)),
  f = exp(bf_1 * a + bf_0),
  u_i = r * a,
  t_i = r + a,
  mat_expr =
    matrix(
      c(
        1 - g_2, 0, v * (1-g_1) * f,
        g_2 * s_1, 0, v * g_1 * s_1 * f,
        0, s_2 * s_3, 0
      ),
      nrow = 3,
      byrow = TRUE

  )
)

# Use make_data_list() for constants and the initial population vector

data <- make_data_list(
  v = 0.8228,
  g_1 = 0.5503,
  g_2 = 0.3171,
  bs2_2 = 0.0016,
  bs2_1 = -0.0664,
  bs2_0 = -0.156,
  bs3_1 = -0.289,
  bf_1 = -0.0389,
  bf_0 = 7.489,
  s_1 = 0.5,
  initial_population_vector = c(s = 10, r = 0, a = 0)
)
```

Note that all constants in the `mat_exprs` list appear in the `data_list`. This includes coefficients in the density dependent expressions (`bs2_0`, `bs2_1`, `bs2_2`, `bs3_0`, `bs3_1` `bf_0`, `bf_1`).

If you are interested in exploring how altering constant vital rates changes the way density dependence plays out, you can loop across the `make_data_list()` call and substitute in new values. You can do the same thing for the coefficients in each density dependent expression.

``` r
source('R/functions.R')

# make this outside the loop, the expressions are not changing
exprs <- make_mat_exprs(
  s_2 = 1/(1 + exp(bs2_2 * u_i + bs2_1 * t_i + bs2_0)),
  s_3 = exp(bs3_1 * log(r + 1)),
  f = exp(bf_1 * a + bf_0),
  u_i = r * a,
  t_i = r + a,
  mat_expr =
    matrix(
      c(
        1 - g_2, 0, v * (1-g_1) * f,
        g_2 * s_1, 0, v * g_1 * s_1 * f,
        0, s_2 * s_3, 0
      ),
      nrow = 3,
      byrow = TRUE

  )
)

# alter s_1 to see how this changes density dependent dynamics
lambdas <- list()
densities <- list()

for(i in seq(0.01, 1, 0.05)) {
  data <- make_data_list(
    v = 0.8228,
    g_1 = 0.5503,
    g_2 = 0.3171,
    bs2_2 = 0.0016,
    bs2_1 = -0.0664,
    bs2_0 = -0.156,
    bs3_1 = -0.289,
    bf_1 = -0.0389,
    bf_0 = 7.489,
    s_1 = i, # note the substitution here
    initial_population_vector = c(s = 10, r = 0, a = 0)
  )
  
  temp_data <- iterate_dd_mat(
    data_list = data,
    mat_exprs = exprs,
    n_generations = 100
  )
  
  lambdas[[paste0('s_1_', i, sep = "")]] <- temp_data$growth_rates$lambda
  densities[[paste0('s_1_', i, sep = "")]] <- apply(temp_data$stage_vectors[ ,2:3],
                                                    1, 
                                                    FUN = sum)
}

# plot the densities. You can change "densities" to "lambdas" in the code
# below to examine growth rates between each iteration
par(mfrow = c(4,5))

lapply(1:20, function(x) {
  plot(densities[[x]], 
       main = names(densities)[x],
       type = 'l')
  
  })
```

### Environmentally-dependent matrices

Iain Stott already has an implementation of these in [popdemo](https://github.com/iainmstott/popdemo). His implementation covers instances where multiple matrices are parameterized with constant values and then are selected either randomly, with Markov Chains, or a user-supplied sequence. I'm working on integrating the density dependent code into the `Projection` class that he wrote for these so that there is a common-ish feel to working with these stochastic matrices.

Still to be implemented are matrices where individual elements are functions of some environmental variable. These will be very similar to the density-dependent matrices, but I'm toying with the idea of adding an additional `env_data` list to the structure. Unfortunately, I have yet to find a paper that reports *all* of the parameters I need to re-create them, so I'm putting this on hold for now to focus on getting the density dependent matrices working across all/as many cases as possible.
