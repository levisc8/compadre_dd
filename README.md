
Attempting to create a framework for density dependence in Comapdre
-------------------------------------------------------------------

This is pretty experimental, but the idea is to store matrix elements (which may be constants *or* expressions) in a long data format, then use a couple functions to take a user- or database supplied population vector and generate a density dependent matrix. Outputs are then generated via iteration.

### Density-dependent matrices

These are now implemented. `iterate_dd_mat()` can handle both user-supplied and data base matrices. Additionally, `make_mat_exprs()` is now smart enough to know when to wrap elements in calls to `eval_tidy()` and `quo()` so that end users and programmers do not need to fully understand how/why these are being used.

`CompadreDDM` matrices do not look like other matrices stored in Compadre. They are lists with 2 elements: a `data_list` which contains values for each parameter and a `mat_exprs` list which contains expressions to calculate density dependent vital rates (e.g. survival, growth, reproduction). Additionally, the `mat_exprs` list contains an expression for the matrix (`mat_expr`). Hopefully, the example below clarifies how these work.

``` r
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

# Use this for constants and the initial population vector
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

If users are interested in exploring how altering constant vital rates alters the way density dependence plays out, this is as simple as looping across the `make_data_list()` call and substituting in new values.

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
    s_1 = i,
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

par(mfrow = c(4,5))

lapply(1:20, function(x) {
  plot(densities[[x]], 
       main = names(densities)[x],
       type = 'l')
  
  })
```
