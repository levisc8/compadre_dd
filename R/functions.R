# Functions to potentially include in Rcompadre to enable working with density
# dependent matrices. Basically, my idea is to modify the relevant mat entries
# so that instead of a CompadreM object in the row, we have a CompadreDDM or something like
# that. It would contain additional information needed to iterate the matrix,
# generate functional forms, etc. The functions depend on rlang which has no
# additional dependencies (except R >= 3.1.0)

# mat_exprs: Must contain all density dependent expressions and an expression
# to generate the matrix itself. terms in density depedendent expressions that
# are themselves functions of something else must be wrapped in calls to
# rlang::eval_tidy and all expressions must be wrapped in calls to rlang::quo.
# Additionally, all density dependent terms must be wrapped in eval_tidy in the
# matrix expression itself.

library(rlang)

iterate_dd_mat <- function(...,
                           n_generations,
                           target_output,
                           init_p_vec) {

  UseMethod('iterate_dd_mat')
}

iterate_dd_mat.CompadreDDM <- function(dd_mat_data,
                                       n_generations,
                                       target_output = c('stage_vectors',
                                                         'growth_rates'),
                                       init_p_vec = NULL) {

  # potential to make dd_mat_data[[fun_arg]] so you can choose which matrix
  # to iterate (U, F, C). For now, just implement for matA

  dd_mat_data <- dd_mat_data$matA

  .iterate_dd_mat_impl(dd_mat_data,
                       n_generations,
                       target_output,
                       init_p_vec)

}

iterate_dd_mat.list <- function(data_list,
                                mat_exprs,
                                n_generations,
                                target_output = c('stage_vectors',
                                                    'growth_rates'), # for now
                                init_p_vec = NULL) {

  dd_mat_data <- list(data_list = data_list,
                      mat_exprs = mat_exprs)

  .iterate_dd_mat_impl(dd_mat_data,
                       n_generations,
                       target_output,
                       init_p_vec)

}

.iterate_dd_mat_impl <- function(mat_data,
                                 n_generations,
                                 target_output,
                                 init_p_vec) {

  data_list <- mat_data$data_list
  mat_exprs <- mat_data$mat_exprs

  if(is.null(init_p_vec)) {
    init_p_vec <- data_list$p_vec
  } else {
    data_list$p_vec <- init_p_vec
  }

  # insulated environment for the iterations to take place
  eval_env <- rlang::env()

  # Set this as the environment for each quosure in mat_exprs
  lapply(mat_exprs, function(x) quo_set_env(x, eval_env))
  # add data
  rlang::env_bind(eval_env,
                  !!! data_list,
                  !!! init_p_vec)

  # create lazy bindings to the evaluation environment
  rlang::env_bind_exprs(eval_env,
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
    new_p_vec <- as.numeric(rlang::eval_tidy(eval_env$mat_expr) %*% eval_env$p_vec)

    # assign stage names to new population vector
    names(new_p_vec) <- names(eval_env$p_vec)

    # bind new names to eval_env for next iteration
    rlang::env_bind(eval_env,
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

    out_obj$growth_rates$lambda[iteration] <-
      sum(to_insert) /
      out_obj$growth_rates$n_tot[(iteration - 1)]
  }

  return(out_obj)

}


# ...: named parameters of the matrix. Every parameter in the functions in
# make_mat_exprs should appear in here. Additionally, constants used in other
# parts of the matrix should also be supplied here
# initial_population_vector: can be user-supplied or pulled from the data base
# (if we decide to store those). otherwise, we will have to decide on sensible
# defaults to use.

make_data_list <- function(...,
                           initial_population_vector = NULL) {
  data <- list(...)
  data$p_vec <- initial_population_vector

  return(data)
}

# These are the functions that correspond to density dependent matrix elements.
# These should contain all of the equations as quosures. Additionally, elements
# that are themselves functions of something else should be wrapped in eval_tidy
# calls
make_mat_exprs <- function(...) {

  mat_exprs <- enquos(...)

  text_list <- wrap_eval_tidys(mat_exprs)

  out <- lapply(text_list, wrap_quos)

  return(out)
}

wrap_eval_tidys <- function(mat_exprs) {

  # turn quos into strings, extract the LHS variable names
  text_list <- lapply(mat_exprs, rlang::quo_text)
  LHS <- names(text_list)

  for(i in seq_along(text_list)) {
    for(j in seq_along(LHS)){

      # substitute eval_tidy(var) if anything in right hand side appears in
      # left hand side
      reg_expression <- make_mat_regexpr(LHS[j])

      text_list[[i]] <- gsub(reg_expression,
                             paste0('eval_tidy(', LHS[j], ')'),
                             text_list[[i]])

    }
  }


  return(text_list)
}

make_mat_regexpr <- function(var) {
    paste0('(\\b', var, '\\b)')
}

wrap_quos <- function(mat_exprs_w_evals) {

  string <- paste0('quo(', mat_exprs_w_evals, ')')
  out <- rlang::parse_expr(string)
  enquo(out)

}
