# BLHS

calc_blhs_energy <- function(mat, lambda) {
  if (nrow(mat) <= 1L || ncol(mat) == 0L)
    return(0)
  dist_mat <- as.matrix(dist(mat, method = "manhattan"))
  lower_idx <- lower.tri(dist_mat)
  distances <- dist_mat[lower_idx]
  distances[distances < 1e-12] <- 1e-12
  dist_sum <- sum((ncol(mat) / distances)^lambda)
  num_pairs <- nrow(mat) * (nrow(mat) - 1L) / 2L
  (dist_sum / num_pairs)^(1 / lambda)
}

calc_blhs_correlation <- function(mat) {
  if (ncol(mat) < 2L)
    return(1e-10)
  corr_mat <- suppressWarnings(cor(mat, use = "pairwise.complete.obs",
    method = "pearson"))
  corr_mat[is.na(corr_mat)] <- 0
  lower_idx <- lower.tri(corr_mat)
  rho2_mean <- mean(corr_mat[lower_idx]^2, na.rm = TRUE)
  if (is.na(rho2_mean) || rho2_mean <= 1e-12) 1e-12 else rho2_mean
}

evaluate_blhs_criteria <- function(mat, lambda) {
  if (ncol(mat) == 1L)
    return(calc_blhs_energy(mat, lambda))
  phi <- calc_blhs_energy(mat, lambda)
  rho2 <- calc_blhs_correlation(mat)
  phi + log10(rho2)
}

scale_matrix <- function(mat) {
  u01 <- (mat - 1 + runif(length(mat))) / nrow(mat)
  2 * u01 - 1
}

optimize_blhs <- function(n, d, lambda, fixed_subspace = NULL, budget = 500L, patience = 50L, tol = 1e-6) {
  if (n <= 1L || d <= 0L)
    return(matrix(runif(n * d), nrow = n, ncol = d))
  if (d == 1L && is.null(fixed_subspace))
    return(matrix((sample.int(n) - 1 + runif(n)) / n, ncol = 1L))
  current_perm <- matrix(0L, nrow = n, ncol = d)
  for (j in seq_len(d))
    current_perm[, j] <- sample.int(n)
  current_scaled <- scale_matrix(current_perm)
  full_current <- if (!is.null(fixed_subspace)) cbind(fixed_subspace, current_scaled) else current_scaled
  current_fit <- evaluate_blhs_criteria(full_current, lambda)
  no_change_count <- 0L
  for (eval in seq_len(budget)) {
    cand_perm <- current_perm
    n_mutate <- max(1L, rbinom(1L, size = d, prob = 1 / d))
    cols_idx <- sample.int(d, size = n_mutate)
    for (col in cols_idx) {
      swap_rows <- sample.int(n, size = 2L)
      tmp <- cand_perm[swap_rows[1L], col]
      cand_perm[swap_rows[1L], col] <- cand_perm[swap_rows[2L], col]
      cand_perm[swap_rows[2L], col] <- tmp
    }
    cand_scaled <- scale_matrix(cand_perm)
    full_cand <- if (!is.null(fixed_subspace)) cbind(fixed_subspace, cand_scaled) else cand_scaled
    cand_fit <- evaluate_blhs_criteria(full_cand, lambda)
    if (cand_fit < (current_fit - tol)) {
      current_perm <- cand_perm
      current_fit <- cand_fit
      current_scaled <- cand_scaled
      no_change_count <- 0L
    } else {
      no_change_count <- no_change_count + 1L
      if (cand_fit <= current_fit) {
        current_perm <- cand_perm
        current_fit <- cand_fit
        current_scaled <- cand_scaled
      }
    }
    if (!is.null(patience) && no_change_count >= patience)
      break
  }
  (current_perm - 1 + runif(length(current_perm))) / n
}

generate_blhs <- function(parameters, nbConfigurations, repair = NULL, budget = 500L) {
  newConfigurations <- configurations_alloc(parameters[["names"]],
    nrow = nbConfigurations, parameters = parameters)
  global_lambda <- sum(sapply(parameters$get_ordered(), function(p) p[["type"]] %in% c("r", "i"))) + 1L
  u_matrix <- matrix(NA_real_, nrow = nbConfigurations, ncol = length(parameters[["names"]]),
    dimnames = list(NULL, parameters[["names"]]))
  cond_groups <- list()
  for (param in parameters$get_ordered()) {
    cond <- param[["condition"]]
    key <- paste(deparse(cond), collapse = " ")
    if (is.null(cond_groups[[key]]))
      cond_groups[[key]] <- list(condition = cond, params = list())
    cond_groups[[key]]$params[[param[["name"]]]] <- param
  }
  assigned_vars <- character(0L)
  for (group in cond_groups) {
    cond <- group$condition
    idx <- which_satisfied(newConfigurations, cond)
    if (length(idx) == 0L)
      next
    params_group <- group$params
    variable_params <- list()
    for (pname in names(params_group)) {
      param <- params_group[[pname]]
      if (param[["isFixed"]]) {
        set(newConfigurations, i = idx, j = pname, value = param[["domain"]])
      } else {
        variable_params[[pname]] <- param
      }
    }
    if (length(variable_params) == 0L)
      next
    fixed_subspace <- if (length(assigned_vars) > 0L) u_matrix[idx, assigned_vars, drop = FALSE] else NULL
    u_block <- optimize_blhs(length(idx), length(variable_params), global_lambda, fixed_subspace, budget)
    u_matrix[idx, names(variable_params)] <- 2 * u_block - 1
    col_idx <- 1L
    for (pname in names(variable_params)) {
      param <- variable_params[[pname]]
      u_col <- u_block[, col_idx]
      if (param[["is_dependent"]]) {
        dep_cols <- parameters$depends[[pname]]
        dep_subset <- as.data.frame(newConfigurations[idx, , drop = FALSE])
        vals <- sapply(seq_len(length(idx)), function(k) {
          dep_dom <- get_dependent_domain(param, dep_subset[k, , drop = FALSE])
          param_quantile(param, probs = u_col[k], domain = dep_dom)
        })
        set(newConfigurations, i = idx, j = pname, value = vals)
      } else {
        set(newConfigurations, i = idx, j = pname, value = param_quantile(param, probs = u_col))
      }
      col_idx <- col_idx + 1L
    }
    assigned_vars <- c(assigned_vars, names(variable_params))
  }
  repairConfigurations(newConfigurations, parameters, repair)
  newConfigurations
}

sampleBLHS <- function(parameters, nbConfigurations, repair = NULL) {
  newConfigurations <- generate_blhs(parameters, nbConfigurations, repair)
  forbidden <- parameters$forbidden
  if (!is.null(forbidden)) {
    retries <- 100L
    repeat {
      newConfigurations <- filter_forbidden(newConfigurations, forbidden)
      needed <- nbConfigurations - nrow(newConfigurations)
      if (needed == 0L)
        break
      newConfigurations <- rbindlist(list(newConfigurations,
        generate_blhs(parameters, needed, repair = repair)))
      retries <- retries - 1L
      if (retries == 0L) {
        irace_error("irace tried 100 times to sample a valid configuration with BLHS without success.")
      }
    }
  }
  set(newConfigurations, j = ".PARENT.", value = NA_integer_)
  newConfigurations
}


# Chaotic with Tent Map

tent_map <- function(n, x0 = NULL, mu = 1.9999) {
  if (n <= 0L) return(numeric(0))
  if (is.null(x0)) {
    x0 <- runif(1L, 0.001, 0.999)
    while (x0 %in% c(0.0, 0.25, 0.5, 0.75, 1.0)) {
      x0 <- runif(1L, 0.001, 0.999)
    }
  }
  x <- numeric(n)
  xk <- x0
  for (k in seq_len(n)) {
    xk <- if (xk < 0.5) mu * xk else mu * (1 - xk)
    if (xk <= 0 || xk >= 1)
      xk <- runif(1L, 0.001, 0.999)
    x[k] <- xk
  }
  x
}

generate_chaotic <- function(parameters, nbConfigurations, repair = NULL) {
  newConfigurations <- configurations_alloc(parameters[["names"]],
    nrow = nbConfigurations, parameters = parameters)
  for (param in parameters$get_ordered()) {
    pname <- param[["name"]]
    idx <- which_satisfied(newConfigurations, param[["condition"]])
    if (length(idx) == 0L)
      next
    if (param[["isFixed"]]) {
      set(newConfigurations, i = idx, j = pname, value = param[["domain"]])
      next
    }
    if (param[["is_dependent"]]) {
      sample_dep_chaos <- function(x)
        param_quantile(param, probs = tent_map(1L), domain = get_dependent_domain(param, x))
      newConfigurations[idx, c(pname) := sample_dep_chaos(.SD), by=.I, .SDcols=parameters$depends[[pname]]]
    } else {
      set(newConfigurations, i = idx, j = pname, value = param_quantile(param, probs = tent_map(length(idx))))
    }
  }
  repairConfigurations(newConfigurations, parameters, repair)
  newConfigurations
}

sampleChaotic <- function(parameters, nbConfigurations, repair = NULL) {
  newConfigurations <- generate_chaotic(parameters, nbConfigurations, repair)
  forbidden <- parameters$forbidden
  if (!is.null(forbidden)) {
    retries <- 100L
    repeat {
      newConfigurations <- filter_forbidden(newConfigurations, forbidden)
      needed <- nbConfigurations - nrow(newConfigurations)
      if (needed == 0L)
        break
      newConfigurations <- rbindlist(list(newConfigurations,
        generate_chaotic(parameters, needed, repair = repair)))
      retries <- retries - 1L
      if (retries == 0L) {
        irace_error("irace tried 100 times to sample chaocticly a configuration not forbidden without success, perhaps your constraints are too strict?")
      }
    }
  }
  set(newConfigurations, j = ".PARENT.", value = NA_integer_)
  newConfigurations
}


# New Sobol Sequence

generate_new_sobol <- function(parameters, nbConfigurations, repair = NULL) {
  newConfigurations <- configurations_alloc(parameters[["names"]],
    nrow = nbConfigurations, parameters = parameters)
  cond_groups <- list()
  for (param in parameters$get_ordered()) {
    cond <- param[["condition"]]
    key <- paste(deparse(cond), collapse = " ")
    if (is.null(cond_groups[[key]]))
      cond_groups[[key]] <- list(condition = cond, params = list())
    cond_groups[[key]]$params[[param[["name"]]]] <- param
  }
  for (group in cond_groups) {
    cond <- group$condition
    idx <- which_satisfied(newConfigurations, cond)
    if (length(idx) == 0L)
      next
    params_group <- group$params
    variable_params <- list()
    for (pname in names(params_group)) {
      param <- params_group[[pname]]
      if (param[["isFixed"]]) {
        set(newConfigurations, i = idx, j = pname, value = param[["domain"]])
      } else {
        variable_params[[pname]] <- param
      }
    }
    if (length(variable_params) == 0L)
      next
    u_sobol <- spacefillr::generate_sobol_set(length(idx), dim = length(variable_params),
      seed = runif_integer(size = 1L))
    col_idx <- 1L
    for (pname in names(variable_params)) {
      param <- variable_params[[pname]]
      u_col <- u_sobol[, col_idx]
      if (param[["is_dependent"]]) {
        dep_cols <- parameters$depends[[pname]]
        dep_subset <- as.data.frame(newConfigurations[idx, , drop = FALSE])
        vals <- sapply(seq_len(length(idx)), function(k) {
          dep_dom <- get_dependent_domain(param, dep_subset[k, , drop = FALSE])
          param_quantile(param, probs = u_col[k], domain = dep_dom)
        })
        set(newConfigurations, i = idx, j = pname, value = vals)
      } else {
        set(newConfigurations, i = idx, j = pname, value = param_quantile(param, probs = u_col))
      }
      col_idx <- col_idx + 1L
    }
  }
  repairConfigurations(newConfigurations, parameters, repair)
  newConfigurations
}

sampleNewSobol <- function(parameters, nbConfigurations, repair = NULL) {
  newConfigurations <- generate_new_sobol(parameters, nbConfigurations, repair)
  newConfigurations <- unique(newConfigurations)
  forbidden <- parameters$forbidden
  newConfigurations <- filter_forbidden(newConfigurations, forbidden)
  have <- nrow(newConfigurations)
  if (have < nbConfigurations) {
    needed <- max(ceiling(nbConfigurations + (nbConfigurations - have) * (2 - have / nbConfigurations)),
      min(parameters$nbVariable * 5L, 100L))
    newConfigurations <- generate_new_sobol(parameters, needed, repair)
    newConfigurations <- unique(newConfigurations)
    newConfigurations <- filter_forbidden(newConfigurations, forbidden)
    if (nrow(newConfigurations) == 0L) {
      irace_error("irace tried to sample a configuration not forbidden without success, perhaps your constraints are too strict?")
    }
    newConfigurations <- truncate_rows(newConfigurations, nbConfigurations)
  }
  set(newConfigurations, j = ".PARENT.", value = NA_integer_)
  newConfigurations
}


# Initialization methods

sample_configurations <- function(parameters, nbConfigurations, method, repair = NULL) {
  allConfigurations <- switch(method,
    "blhs" = sampleBLHS(parameters, nbConfigurations, repair = repair),
    "chaotic" = sampleChaotic(parameters, nbConfigurations, repair = repair),
    "sobol" = sampleNewSobol(parameters, nbConfigurations, repair = repair),
    "uniform" = sampleUniform(parameters, nbConfigurations, repair = repair)
  )
  allConfigurations
}