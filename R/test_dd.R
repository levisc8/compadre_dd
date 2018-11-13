# Test Compadre DD ideas

library(rlang)


# Now, let's try w/ density dependent expressions
# Now, for the density dependent expressions:
# s_2 = survival over the summer of rosettes
# s_3 = survival over the winter of rosettes
# u = interaction between rosettes and adults (a * r) (from population vector)
# a = density of adults (from population vector)
# r = density of rosettes (from population vector)
# b* = regression coefficients for each term
# s_1 = a range between 0.1 and 0.9. Apparently they couldn't test this

s_2 <- quo(1/(1 + exp(bs2_2 * eval_tidy(u_i) + bs2_1 * eval_tidy(t_i) + bs2_0)))

s_3 <- quo(exp(bs3_1 * log(r + 1)))

f <- quo(exp(bf_1 * a + bf_0))

u_i <- quo(r * a)
t_i <- quo(r + a)

v <- 0.8228
g_1 <- 0.5503
g_2 <- 0.3171
bs2_2 <- 0.0016
bs2_1 <- -0.0664
bs2_0 <- -0.156
bs3_1 <- -0.289
bf_1 <- -0.0389
bf_0 <- 7.489
s_1 <- 0.5

data_mat <- quo(
  matrix(
    c(
      1 - g_2, 0, v * (1-g_1) * eval_tidy(f),
      g_2 * s_1, 0, v * g_1 * s_1 * eval_tidy(f),
      0, eval_tidy(s_2) * eval_tidy(s_3), 0
    ),
    nrow = 3,
    byrow = TRUE
  )
)

a <- 0
r <- 0
s <- 10

p_vec <- c(s,r,a)

eval_tidy(f) # make sure it works
n_p_vec <- eval_tidy(data_mat) %*% p_vec # iterate the matrix once

# iterate it again
a <- n_p_vec[3]
r <- n_p_vec[2]
s <- n_p_vec[1]

n_p_vec_2 <- eval_tidy(data_mat) %*% n_p_vec

# doing this by hand is lame, lets loop it for 100 generations
n_gen <- 100

output <- data.frame(s = c(10, rep(NA_real_, n_gen)),
                     r = c(0, rep(NA_real_, n_gen)),
                     a = c(0, rep(NA_real_, n_gen)),
                     n_tot = c(10, rep(NA_real_, n_gen)),
                     g_rate = rep(NA_real_, n_gen + 1))

for(i in seq_len(n_gen)) {
  p_vec <- as.numeric(output[i, 1:3])
  new_p_vec <- eval_tidy(data_mat) %*% p_vec
  output[(i+1), 1:3] <- new_p_vec
  a <- new_p_vec[3]
  r <- new_p_vec[2]
  s <- new_p_vec[1]
  output[(i+1), 4] <- sum(new_p_vec)
  output[(i+1), 5] <- sum(new_p_vec)/output[i, 4]

}

par(mfrow = c(2,3))
for(i in seq_along(output)) {
  plot(output[ ,i])
}
