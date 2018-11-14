# Functions to potentially include in Rcompadre to enable working with density
# dependent matrices. Basically, my idea is to restructure the mat entries
# so that instead of a list of CompadreM, we have CompadreDDM or something like
# that. It would contain additional information needed to iterate the matrix,
# generate functional forms, etc. The functions depend on rlang which has no
# additional dependencies (except R >= 3.1.0)

# mat_exprs: Must contain all density dependent expressions and an expression
# to generate the matrix itself. terms in density depedendent expressions that
# are themselves functions of something else must be wrapped in calls to
# rlang::eval_tidy and all who expressions must be wrapped in calls to rlang::quo

library(rlang)

iterate_dd_mat <- function(mat_exprs,
                           data_list,
                           n_generations,
                           target_output = c('stage_vectors',
                                             'growth_rates'), # for now
                           init_p_vec = NULL) {

  if(is.null(init_p_vec)) {
    init_p_vec <- data_list$p_vec
  }
  # insulated environment for the iterations to take place
  eval_env <- env()
  # add data and expressions
  env_bind(eval_env,
           !!! data_list,
           !!! data_list$p_vec)

  env_bind_exprs(eval_env,
                 !!! mat_exprs,
                 .eval_env = eval_env)

  # set up place to store outputs of interest
  output <- list()
  if('stage_vectors' %in% target_output) {
    output$stage_vectors <- data.frame(t(init_p_vec),
                                       n_tot = sum(init_p_vec))
    output$stage_vectors[2:n_generations, ] <- NA_real_
  }

  if('growth_rates' %in% target_output) {
    output$growth_rates <- data.frame(n_tot = c(sum(init_p_vec),
                                                    rep(NA_real_, n_generations - 1)),
                                      lambda = rep(NA_real_, n_generations))
  }

  # iterate!
  for(i in seq(2, n_generations, 1)) {
    # iterate the matrix and generate new population vector
    new_p_vec <- as.numeric(eval_tidy(eval_env$mat_expr) %*% eval_env$p_vec)

    # assign stage names to new population vector
    names(new_p_vec) <- names(eval_env$p_vec)

    # bind new names to eval_env for next iteration
    env_bind(eval_env,
             !!! new_p_vec)
    eval_env$p_vec <- new_p_vec

    # stash outputs in our list
    output <- update_dd_outputs(output, new_p_vec, i, target_output)

  }

  return(output)
}

update_dd_outputs <- function(out_obj,
                              to_insert,
                              iteration,
                              targets) {

  if('stage_vectors' %in% targets) {
    out_obj$stage_vectors[(iteration), 1:(length(to_insert))] <- to_insert
    out_obj$stage_vectors$n_tot[(iteration)] <- sum(to_insert)
  }

  if('growth_rates' %in% targets) {
    out_obj$growth_rates$n_tot[iteration] <- sum(to_insert)

    out_obj$growth_rates$lambda[iteration] <- sum(to_insert) /
                                              out_obj$growth_rates$n_tot[(iteration - 1)]
  }

  return(out_obj)

}


# ...: named constant parameters of the matrix. Generally, I think it should follow the
# convention of a11, a12, a21, a22, etc, but we can consider other options.
# initial_population_vector: can be user-supplied or pulled from the data base
# (if we decide to store those). otherwise, we will have to decide on sensible
# defaults to use.

make_data_list <- function(...,
                           initial_population_vector = NULL) {
  data <- list(...)
  data$p_vec <- initial_population_vector

  return(data)
}

make_mat_exprs <- function(...) {

  out <- enquos(...)

  return(out)
}
