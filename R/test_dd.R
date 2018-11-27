# Test Compadre DD ideas

# rm(list = ls())
library(rlang)


# Let's try w/ density dependent expressions. this is based on Pardini et al 2009
# Complex population dynamics and control of the invasive biennial Alliaria petiolata,
# Ecological Applications

source('R/functions.R')

# Here are the expressions. These would need to be stored in Compadre somewhere,
# as would the data list
# Now, for the density dependent expressions:
# s_2 = survival over the summer of rosettes
# s_3 = survival over the winter of rosettes
# u = interaction between rosettes and adults (a * r) (from population vector)
# a = density of adults (from population vector)
# r = density of rosettes (from population vector)
# b* = regression coefficients for each term
# s_1 = a range between 0.1 and 0.9. Apparently they couldn't test this

# mat_exprs holds the density dependent expressions and the actualy matrix itself.
# make_mat_exprs is now smart enough figure out who needs quoting and evaluating calls,
# so Compadrinos will not need to understand tidy_eval provided the expressions
# and data are correctly specified

# exprs <- make_mat_exprs(
#   s_2 = quo(1/(1 + exp(bs2_2 * eval_tidy(u_i) + bs2_1 * eval_tidy(t_i) + bs2_0))),
#   s_3 = quo(exp(bs3_1 * log(r + 1))),
#   f = quo(exp(bf_1 * a + bf_0)),
#   u_i = quo(r * a),
#   t_i = quo(r + a),
#   mat_expr = quo(
#     matrix(
#       c(
#         1 - g_2, 0, v * (1-g_1) * eval_tidy(f),
#         g_2 * s_1, 0, v * g_1 * s_1 * eval_tidy(f),
#         0, eval_tidy(s_2) * eval_tidy(s_3), 0
#       ),
#       nrow = 3,
#       byrow = TRUE
#     )
#   )
# )

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

# load in old compadre data set
library(RcompadreTidy)
data("Compadre")

dd_mat <- list(matA = list(mat_exprs = exprs,
                           data_list = data))

class(dd_mat) <- 'CompadreDDM'

Compadre@data[151, ] <- NA
Compadre@data$matdata[151] <- dd_mat

hold_breath <- iterate_dd_mat(dd_mat,
                              n_generations = 100)

# Some weird stuff going on w/ nesting it in Compadre. Obviously, this needs
# to be solved.
# test <- iterate_dd_mat(Compadre@data$matdata[[151]],
#                        n_generations = 100)

complete_data <- do.call(cbind, hold_breath)
par(mfrow = c(3,2))

for(i in seq_along(complete_data)) {
  plot(complete_data[ ,i], type = 'l', main = names(complete_data)[i])
}

# Now, try it with user defined list (the same as above, but a different interface)

hold_breath_again <- iterate_dd_mat(data_list = data,
                                    mat_exprs = exprs,
                                    n_generations = 100)

# moment of truth!
identical(hold_breath, hold_breath_again)


