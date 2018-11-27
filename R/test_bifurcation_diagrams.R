# Bifurcation diagrams in appendix B of pardini et al 2009

source('R/functions.R')
library(dplyr)
library(tidyr)
library(ggplot2)

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

# reformat data. I'm taking generations 50 - 100 for this, but I think
# Pardini et al only use generation 50 for their bifurcation plots? It's kind
# of unclear from the text.

for_plot <- densities %>%
  lapply(FUN = function(x) x[50:100]) %>%
  as_tibble() %>%
  gather() %>%
  mutate(s_1 = stringr::str_remove(.$key, 's_1_') %>% as.numeric)

ggplot(for_plot, aes(x = s_1, y = value)) +
  geom_point()

# looks like we get chaos where they get chaos, and the stable cycles
# are almost exactly the same as theirs
