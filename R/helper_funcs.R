# Main HierNet test statistic as implemented in Equation 5 and 11
# Additional inputs:
# Analysis: if TRUE function also returns which interactions are strongest
# Forced: if non-null takes input of which indexes of X matrix we wish to "force" as main effect (see Section 5.2 for further details)
# Group: list (total length of unlisted list should match idx) that "groups" up all relevant interested indexes in idx to implement Equation 5 (only necessary when idx contains indexes we do not want to average the respective effects for and compare together as done in Equation 5)
# Immigration Example: idx = indexes of matrix X that refer to Mexico or European (length of 2)
# Gender Example: idx = indexes of matrix X that refers to Male or Female and all forced Male/Female interactions with Party
# Gender Example: forced = indexes of matrix X of all forced interactions with party affiliation and party affiliation main effects
# Gender Example: group = list(c(1,2), c(3, 4), c(5, 6), c(7,8), c(9, 10)), where c(1,2) for example compares all main and interaction of just male and female. c(3,4) compares all main and interaction effects of male and one factor level from party affiliation, etc.
hiernet_group = function(hiernet_object, idx, X, analysis = 0, group) {
  main = hiernet_object$bp[idx] - hiernet_object$bn[idx]
  main_contribution = vector()
  for (i in 1:length(group)) {
    main_contribution[i] = sum((main[group[[i]]] - mean(main[group[[i]]]))^2)
  }

  I_int = list()
  for (i in 1:length(idx)) {
    I_int[[i]] = (hiernet_object$th[idx[i], ] + hiernet_object$th[, idx[i]])/2
  }


  int_cont_bygroup = vector()
  tracking_all_ints = list()
  for (i in 1:length(group)) {
    in_int = I_int[group[[i]]]
    int_means = Reduce("+", in_int)/length(in_int)
    getting_int_cont = vector()
    for (j in 1:length(in_int[[1]])) {
      rel_ints = sapply(in_int, "[[", j)
      getting_int_cont[j] = sum((rel_ints - int_means[j])^2)
    }
    tracking_all_ints[[i]] = getting_int_cont
    int_cont_bygroup[i] = sum(getting_int_cont)
  }

  # divide by 2 to account for counting both left and right profile
  ts = sum(main_contribution)/2 + sum(int_cont_bygroup)/2

  if (analysis > 0) {
    main_contribution = sum(main_contribution)/2
    int_contribution = sum(int_cont_bygroup)/2
    largest_ints = largest_ints_idx = list()
    for (i in 1:length(group)) {
      largest_ints[[i]] = sort(tracking_all_ints[[i]], decreasing = TRUE)[c(1:(2*analysis))]
      rel_int_name = paste0(colnames(X)[idx][group[[i]]], collapse = ", ")
      largest_ints_idx[[i]] = paste0("interactions with ",colnames(X)[order(tracking_all_ints[[i]], decreasing = TRUE)[c(1:(2*analysis))]])

    }

    largest_ints = unlist(largest_ints)
    largest_ints_idx = unlist(largest_ints_idx)

    #take first and third index (since the second one is a repeat of the first)
    largest_ints_final = sort(largest_ints, decreasing = TRUE)[seq(1, (2*analysis), by = 2)]

    largest_int_contributer = largest_ints_idx[order(largest_ints, decreasing = TRUE)[seq(1, (2*analysis), by = 2)]]

    analysis_list = list()
    analysis_list$observed_TS = ts
    analysis_list$main_contirbution = main_contribution
    analysis_list$int_contribution = int_contribution
    analysis_list$largest_int_size = largest_ints_final
    analysis_list$largest_int_contributer = largest_int_contributer


    return(analysis_list)

  } else {
    return(ts)
  }
}

# Implementing Test statistic to test no profile order effect (equation in Section 3.5)
PO_stat = function(hiernet_object, in_idx_left, in_idx_right, in_idx_respondent) {
  main_1 = hiernet_object$bp[in_idx_left] - hiernet_object$bn[in_idx_left]
  main_2 = hiernet_object$bp[in_idx_right] - hiernet_object$bn[in_idx_right]

  # interaction effects
  within_int_list_left = list()
  within_int_list_right = list()

  for (i in 1:length(in_idx_left)) {
    within_int_list_left[[i]] = (hiernet_object$th[in_idx_left[i], in_idx_left] + hiernet_object$th[in_idx_left, in_idx_left[i]])/2
    within_int_list_right[[i]] = (hiernet_object$th[in_idx_right[i], in_idx_right] + hiernet_object$th[in_idx_right, in_idx_right[i]])/2
  }

  within_diff = unlist(Map("+", within_int_list_left, within_int_list_right))

  #between profiles
  between_int_list_left = list()
  between_int_list_right = list()

  for (i in 1:length(in_idx_left)) {
    between_int_list_left[[i]] = (hiernet_object$th[in_idx_left[i], in_idx_right] + hiernet_object$th[in_idx_right, in_idx_left[i]])/2
    between_int_list_right[[i]] = (hiernet_object$th[in_idx_right[i], in_idx_left] + hiernet_object$th[in_idx_left, in_idx_right[i]])/2
  }

  between_diff = unlist(Map("+", between_int_list_left, between_int_list_right))


  #respondent
  R_int_list_left = list()
  R_int_list_right = list()

  for (i in 1:length(in_idx_left)) {
    R_int_list_left[[i]] = (hiernet_object$th[in_idx_left[i], in_idx_respondent] + hiernet_object$th[in_idx_respondent, in_idx_left[i]])/2
    R_int_list_right[[i]] = (hiernet_object$th[in_idx_right[i], in_idx_respondent] + hiernet_object$th[in_idx_respondent, in_idx_right[i]])/2
  }

  respondent_effects = unlist(Map("+", R_int_list_left, R_int_list_right))

  #division of two in the interactions because I overcount
  stat= sum((main_1 + main_2)^2) + sum(within_diff^2) + sum(between_diff^2) + sum((respondent_effects)^2)
  return(stat)
}

# Implementing Test statistic to test carryover effect (equation in Section 3.5)
CO_stat = function(hiernet_object, idx) {
  I_list = list()
  for (i in 1:length(idx)) {
    I_list[[i]] = (hiernet_object$th[idx[i], -idx] + hiernet_object$th[-idx, idx[i]])/2
  }
  ts = sum(unlist(I_list)^2)
  return(ts)
}

get_CRT_pval = function(x, y, xcols, left_idx, right_idx, design, B, num_cores, profileorder_constraint, lambda, non_factor_idx, in_levs, analysis, p, resample_func_1, resample_func_2, tol, resample_X, full_X, restricted_X, left_allowed, right_allowed, forced, speedup, seed, supplyown_resamples, parallel, nfolds) {
  num_x_levs = levels(x[, xcols[1]])

  if (length(unique(sapply(x, class)[!(1:ncol(x)) %in% non_factor_idx])) > 1) stop("factors provided in formula are not all factors please supply non_factor_idx")

  if (!all(in_levs %in% num_x_levs)) stop("in_levs supplied are not levels in X")
  if (length(left_idx) != length(right_idx)) stop("length left idx does not match right idx")

  if (length(xcols) != 2) stop("supply left and right column for factor X")
  if (!is.null(p)) {
    if (!is.null(in_levs)) {
      if (length(p) != length(xcols)) stop("length of supplied probability should match length of interested levels in_levs")

    } else {
      if (length(num_x_levs) != p) stop("length of supplied probability should match length of number of levels in X")

    }
    if (length(p) != length(xcols)) stop("length of supplied probability should match length of interested levels in_levs")
  }


  # forcing no profile order effect constraint
  if (profileorder_constraint) {
    x_df = x
    n = nrow(x_df)
    empty_df = x_df
    y_new = 1 - y
    for (i in 1:length(left_idx)) {
      empty_df[, left_idx[i]] = x_df[, right_idx[i]]
      empty_df[, right_idx[i]] = x_df[, left_idx[i]]
    }

    final_df = rbind(x_df, empty_df)
    final_df$Y = c(y, y_new)

  } else {
    final_df = x
    final_df$Y = y
  }




  if (is.null(in_levs)) {
    X_names = list()
    for (j in 1:length(xcols)) {
      start_name = colnames(x)[xcols[j]]
      resulting_X_name = vector()
      for (i in 1:length(num_x_levs)) {
        resulting_X_name[i] = paste0(start_name, num_x_levs[i])
      }
      X_names[[j]] = resulting_X_name
    }
    all_X_names = unlist(X_names)


  } else {
    #### focus only on the relevant levels given
    X_names = list()
    for (j in 1:length(xcols)) {
      start_name = colnames(x)[xcols[j]]
      resulting_X_name = vector()
      for (i in 1:length(in_levs)) {
        resulting_X_name[i] = paste0(start_name, in_levs[i])
      }
      X_names[[j]] = resulting_X_name
    }
    all_X_names = unlist(X_names)
  }

  # create X matrix and track index for the relevant effects
  if (is.null(forced)) {

    form = formula(paste0("Y ~ . "))
    X = model.matrix(form, final_df, contrasts.arg = lapply(final_df[, -c(non_factor_idx, ncol(final_df))], contrasts, contrasts = FALSE))[, -1]
    idx_in = (1:ncol(X))[colnames(X) %in% all_X_names]
    all_X_names = X_names
  } else {
    x_names = colnames(x)[xcols]
    forced_names = colnames(x)[forced]
    force_syntax = ""
    for (i in 1:length(forced_names)) {
      force_syntax = paste0(force_syntax, " + ", x_names[1], "*", forced_names[i], " + ", x_names[2], "*", forced_names[i])

    }
    force_syntax = substr(force_syntax, 4, nchar(force_syntax))
    form = formula(paste0("Y ~ . + ", force_syntax))
    X = model.matrix(form, final_df, contrasts.arg = lapply(final_df[, -c(non_factor_idx, ncol(final_df))], contrasts, contrasts = FALSE))[, -1]

    X_forced_names = vector()
    forced_levs = unique(x[, forced[1]])

    all_forced_names = vector()
    for (j in 1:length(forced)) {
      start_name = colnames(x)[forced[j]]
      resulting_forced_name = vector()
      for (i in 1:length(forced_levs)) {
        resulting_forced_name[i] = paste0(start_name, forced_levs[i])
      }
      all_forced_names = c(all_forced_names, resulting_forced_name)
    }

    #repeat this for X:Z and Z:X since model.matrix arbitarily assigns them
    all_possible_ints_1 = list()
    for (j in 1:length(all_forced_names)) {
      temp_vec = list()
      for (i in 1:length(X_names)) {
        temp_vec[[i]] = paste0(all_forced_names[j], ":", X_names[[i]])
      }
      all_possible_ints_1 = c(all_possible_ints_1, temp_vec)
    }

    all_possible_ints_2 = list()
    for (j in 1:length(all_forced_names)) {
      temp_vec = list()
      for (i in 1:length(X_names)) {
        temp_vec[[i]] = paste0(X_names[[i]], ":", all_forced_names[j])
      }
      all_possible_ints_2 = c(all_possible_ints_2, temp_vec)
    }

    all_possible_ints = unlist(c(all_possible_ints_1, all_possible_ints_2))

    idx_in = (1:ncol(X))[colnames(X) %in% c(all_X_names, all_possible_ints)]
    all_X_names = c(X_names, all_possible_ints_1, all_possible_ints_2)


  }

  # this is a quick harmless fix to avoid columns that have all zero
  check_cols = sapply(data.frame(X), function(x) length(unique(x)))
  trouble_cols = as.vector(which(check_cols == 1))
  if (length(trouble_cols) > 0) {
    X[, trouble_cols][1, ] = X[, trouble_cols][1, ] + 1e-5
  }

  in_col_names = colnames(X)[idx_in]

  group_idx_track = vector()
  for (i in 1:length(in_col_names)) {
    group_idx_track[i] = which((sapply(sapply(all_X_names, function(x) which(in_col_names[i] == x)), length)) >0)
  }

  group = list()

  unique_group_idx = unique(group_idx_track)
  for (i in 1:length(unique_group_idx)) {
    group[[i]] = which(group_idx_track == unique_group_idx[i])
  }

  # if in_levs is non-null find idx to resample
  if (is.null(in_levs)) {
    idx_1 = idx_2 = 1:nrow(x)
  } else {
    idx_1 = (1:nrow(x))[x[, xcols[1]] %in% c(in_levs)]
    idx_2 = (1:nrow(x))[x[, xcols[2]] %in% c(in_levs)]
  }
  num_sample_1 = length(idx_1)
  num_sample_2 = length(idx_2)

  set.seed(seed)

  # Creates only Z matrix to perform CV on Z
  if (speedup) {
    Z_df = final_df[, -xcols]
    if (is.null(non_factor_idx)) {
      Z_long = model.matrix(Y ~ . , Z_df, contrasts.arg = lapply(Z_df[, 1:(ncol(Z_df) - 1)], contrasts, contrasts = FALSE))[, -1]
    } else {
      non_factor_names = colnames(x)[non_factor_idx]
      Z_nonfactor_idx = (1:ncol(Z_df))[(colnames(Z_df) %in% non_factor_names)]
      Z_long = model.matrix(Y ~ . , Z_df, contrasts.arg = lapply(Z_df[, -c(Z_nonfactor_idx, ncol(Z_df))], contrasts, contrasts = FALSE))[, -1]

    }

    # obtains CV lambda
    best_lam = hierNet_logistic_CV(lambda, nfolds = nfolds, X = Z_long, y_var = Z_df$Y, tol = tol, constraint = profileorder_constraint, seed = seed)

    print("Initial step completed: finished computing cross validated lambda")

  } else {
    best_lam = hierNet_logistic_CV(lambda, nfolds = nfolds, X = as.matrix(X), y_var = final_df$Y, tol = tol, constraint = profileorder_constraint, seed = seed)

  }

  if (speedup) {
    chosen_initial = sample(c(0:B), size = 1)
    if (chosen_initial != 0) {
      x_df = x
      if (is.null(supplyown_resamples)) {
        if (design == "Uniform") {
          newx_1 = factor(resample_func_1(resample_X, nrow(x)))
          newx_2 = factor(resample_func_2(resample_X, nrow(x)))
        }
        if (design == "Constrained Uniform") {
          newx_1 = factor(resample_func_1(full_X, restricted_X, nrow(x), left_allowed))
          newx_2 = factor(resample_func_2(full_X, restricted_X, nrow(x), left_allowed))
        }
        if (design == "Nonuniform") {
          newx_1 = factor(resample_func_1(resample_X, nrow(x), p))
          newx_2 = factor(resample_func_2(resample_X, nrow(x), p))
        }

      } else {
        newx_1 = factor(supplyown_resamples[[chosen_initial]][, 1])
        newx_2 = factor(supplyown_resamples[[chosen_initial]][, 2])
      }

      x_df[, xcols[1]][idx_1] = newx_1[idx_1]
      x_df[, xcols[2]][idx_2] = newx_2[idx_2]
      # forcing no profile order effect constraint
      if (profileorder_constraint) {
        n = nrow(x_df)
        empty_df = x_df
        y_new = 1 - y
        for (i in 1:length(left_idx)) {
          empty_df[, left_idx[i]] = x_df[, right_idx[i]]
          empty_df[, right_idx[i]] = x_df[, left_idx[i]]
        }

        final_df = rbind(x_df, empty_df)
        final_df$Y = c(y, y_new)

      } else {
        final_df = x_df
        final_df$Y = y
      }

      X_initial = model.matrix(form, final_df, contrasts.arg = lapply(final_df[, -c(non_factor_idx, ncol(final_df))], contrasts, contrasts = FALSE))[, -1]
      # this is a quick harmless fix to avoid columns that have all zero
      check_cols = sapply(data.frame(X_initial), function(x) length(unique(x)))
      trouble_cols = as.vector(which(check_cols == 1))
      if (length(trouble_cols) > 0) {
        X_initial[, trouble_cols][1, ] = X_initial[, trouble_cols][1, ] + 1e-5
      }

      invisible(capture.output(initial <- hierNet_logistic(as.matrix(X_initial), final_df$Y, lam= best_lam, tol = tol)))
      ## new changed way
    } else {
      invisible(capture.output(initial <- hierNet_logistic(as.matrix(X), final_df$Y, lam= best_lam, tol = tol)))
    }
    aa = initial
  } else {
    aa = NULL
  }


  invisible(capture.output(fit <- hierNet_logistic(as.matrix(X), final_df$Y, lam= best_lam, tol = tol, aa = aa)))

  obs_test_stat = hiernet_group(fit, idx = idx_in, X = X, analysis = analysis, group = group)

  print("Initial step completed: finished computing observed test statistic. Now entering resampling step.....")

  ## parallel computing setup
  if (parallel) {
    cl <- snow::makeCluster(num_cores)
    doSNOW::registerDoSNOW(cl)
    iterations <- B
    pb <- utils::txtProgressBar(max = iterations, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    ##

    e <- foreach::foreach(j = 1:iterations, .combine = rbind,
                 .options.snow = opts) %dopar%
      {
        set.seed(j + B)
        library(CRTConjoint)

        x_df = x

        if (is.null(supplyown_resamples)) {
          if (design == "Uniform") {
            newx_1 = factor(resample_func_1(resample_X, nrow(x)))
            newx_2 = factor(resample_func_2(resample_X, nrow(x)))
          }
          if (design == "Constrained Uniform") {
            newx_1 = factor(resample_func_1(full_X, restricted_X, nrow(x), left_allowed))
            newx_2 = factor(resample_func_2(full_X, restricted_X, nrow(x), right_allowed))
          }
          if (design == "Nonuniform") {
            newx_1 = factor(resample_func_1(resample_X, nrow(x), p))
            newx_2 = factor(resample_func_2(resample_X, nrow(x), p))
          }

        } else {
          newx_1 = factor(supplyown_resamples[[j]][, 1])
          newx_2 = factor(supplyown_resamples[[j]][, 2])
        }

        x_df[, xcols[1]][idx_1] = newx_1[idx_1]
        x_df[, xcols[2]][idx_2] = newx_2[idx_2]
        # forcing no profile order effect constraint
        if (profileorder_constraint) {
          n = nrow(x_df)
          empty_df = x_df
          y_new = 1 - y
          for (i in 1:length(left_idx)) {
            empty_df[, left_idx[i]] = x_df[, right_idx[i]]
            empty_df[, right_idx[i]] = x_df[, left_idx[i]]
          }

          final_df = rbind(x_df, empty_df)
          final_df$Y = c(y, y_new)

        } else {
          final_df = x_df
          final_df$Y = y
        }

        X = model.matrix(form, final_df, contrasts.arg = lapply(final_df[, -c(non_factor_idx, ncol(final_df))], contrasts, contrasts = FALSE))[, -1]

        # this is a quick harmless fix to avoid columns that have all zero
        check_cols = sapply(data.frame(X), function(x) length(unique(x)))
        trouble_cols = as.vector(which(check_cols == 1))
        if (length(trouble_cols) > 0) {
          X[, trouble_cols][1, ] = X[, trouble_cols][1, ] + 1e-5
        }

        if (!speedup) {
          best_lam = hierNet_logistic_CV(lambda, nfolds = nfolds, X = as.matrix(X), y_var = final_df$Y, tol = tol, constraint = profileorder_constraint, seed = seed)
        }

        invisible(capture.output(fit <- hierNet_logistic(as.matrix(X), final_df$Y, lam= best_lam, tol = tol, aa = aa)))

        resamp_TS = hiernet_group(fit, idx = idx_in, X = X, analysis = 0, group = group)

        return(resamp_TS)
      }

    snow::stopCluster(cl)

  } else {
    e = vector()
    for (j in 1:B) {
      set.seed(j + B)

      x_df = x

      if (is.null(supplyown_resamples)) {
        if (design == "Uniform") {
          newx_1 = factor(resample_func_1(resample_X, nrow(x)))
          newx_2 = factor(resample_func_2(resample_X, nrow(x)))
        }
        if (design == "Constrained Uniform") {
          newx_1 = factor(resample_func_1(full_X, restricted_X, nrow(x), left_allowed))
          newx_2 = factor(resample_func_2(full_X, restricted_X, nrow(x), right_allowed))
        }
        if (design == "Nonuniform") {
          newx_1 = factor(resample_func_1(resample_X, nrow(x), p))
          newx_2 = factor(resample_func_2(resample_X, nrow(x), p))
        }

      } else {
        newx_1 = factor(supplyown_resamples[[j]][, 1])
        newx_2 = factor(supplyown_resamples[[j]][, 2])
      }

      x_df[, xcols[1]][idx_1] = newx_1[idx_1]
      x_df[, xcols[2]][idx_2] = newx_2[idx_2]
      # forcing no profile order effect constraint
      if (profileorder_constraint) {
        n = nrow(x_df)
        empty_df = x_df
        y_new = 1 - y
        for (i in 1:length(left_idx)) {
          empty_df[, left_idx[i]] = x_df[, right_idx[i]]
          empty_df[, right_idx[i]] = x_df[, left_idx[i]]
        }

        final_df = rbind(x_df, empty_df)
        final_df$Y = c(y, y_new)

      } else {
        final_df = x_df
        final_df$Y = y
      }

      X = model.matrix(form, final_df, contrasts.arg = lapply(final_df[, -c(non_factor_idx, ncol(final_df))], contrasts, contrasts = FALSE))[, -1]
      # this is a quick harmless fix to avoid columns that have all zero
      check_cols = sapply(data.frame(X), function(x) length(unique(x)))
      trouble_cols = as.vector(which(check_cols == 1))
      if (length(trouble_cols) > 0) {
        X[, trouble_cols][1, ] = X[, trouble_cols][1, ] + 1e-5
      }

      if (!speedup) {
        best_lam = hierNet_logistic_CV(lambda, nfolds = nfolds, X = as.matrix(X), y_var = final_df$Y, tol = tol, constraint = profileorder_constraint, seed = seed)
      }

      invisible(capture.output(fit <- hierNet_logistic(as.matrix(X), final_df$Y, lam= best_lam, tol = tol, aa = aa)))

      e[j] = hiernet_group(fit, idx = idx_in, X = X, analysis = 0, group = group)
      print(paste0("Done with task: ",j, " out of ", B, " resamples"))
    }


  }
  p_val = (length(which(e >= as.numeric(obs_test_stat[1]))) + 1)/(B + 1)

  out = list()
  out$p_val = p_val
  out$obs_test_stat = obs_test_stat
  out$resampled_test_stat = e
  out$tol = tol
  if (speedup) {
    out$lam = best_lam
  }
  out$hiernet_fit = fit
  out$seed = seed

  return(out)
}

get_profileordereffect = function(x, y, left_idx, right_idx, B, num_cores, lambda, non_factor_idx, tol, speedup, seed, parallel, nfolds) {

  if (unique(sapply(x[, left_idx[!(left_idx %in% non_factor_idx)]], class)) != "factor") stop("left and right factors are not all factors please supply non_factor_idx")


  if (length(left_idx) != length(right_idx)) stop("length left idx does not match right idx")


  final_df = x
  final_df$Y = y

  form = formula(paste0("Y ~ . "))
  X = model.matrix(form, final_df, contrasts.arg = lapply(final_df[, -c(non_factor_idx, ncol(final_df))], contrasts, contrasts = FALSE))[, -1]

  set.seed(seed)
  if (speedup) {
    chosen_initial = sample(c(0:B), size = 1)
    if (chosen_initial != 0) {
      x_df = x

      resampling_df = final_df
      a = sample(c(0,1), size =nrow(x), replace = TRUE)
      sample_idx = which(a == 1)

      kept = resampling_df[-sample_idx, ]

      b = resampling_df[sample_idx, ]
      new_y = 1 - b$Y
      new_first = b[, left_idx]
      new_second = b[, right_idx]
      name_1 = colnames(new_first)
      name_2 = colnames(new_second)
      colnames(new_first) = name_2
      colnames(new_second) = name_1
      new_df = cbind(new_second, new_first)
      new_df = cbind(new_df, b[, -c(left_idx, right_idx, ncol(b))])

      new_df$Y = new_y
      final_df_initial = rbind(kept, new_df)




      X_initial = model.matrix(form, final_df_initial, contrasts.arg = lapply(final_df_initial[, -c(non_factor_idx, ncol(final_df_initial))], contrasts, contrasts = FALSE))[, -1]

      best_lam = hierNet_logistic_CV(lambda, nfolds = nfolds, X = X_initial, y_var = final_df_initial$Y, tol = tol, constraint = FALSE, seed = seed)

      invisible(capture.output(initial <- hierNet_logistic(as.matrix(X_initial), final_df$Y, lam= best_lam, tol = tol)))

      ## new changed way
    } else {

      best_lam = hierNet_logistic_CV(lambda, nfolds = nfolds, X = X, y_var = final_df$Y, tol = tol, constraint = FALSE, seed = seed)

      invisible(capture.output(initial <- hierNet_logistic(as.matrix(X), final_df$Y, lam= best_lam, tol = tol)))
    }
    aa = initial
  } else {

    best_lam = hierNet_logistic_CV(lambda, nfolds = nfolds, X = X, y_var = final_df$Y, tol = tol, constraint = FALSE, seed = seed)

    aa = NULL
  }


  invisible(capture.output(fit <- hierNet_logistic(as.matrix(X), final_df$Y, lam= best_lam, tol = tol, aa = aa)))

  # get relevant index
  left_X = final_df[, c(left_idx, ncol(final_df))]
  X_left = model.matrix(form, left_X, contrasts.arg = lapply(left_X[, -c(non_factor_idx, ncol(left_X))], contrasts, contrasts = FALSE))[, -1]
  in_idx_left =  (1:ncol(X))[colnames(X) %in% colnames(X_left)]

  right_X = final_df[, c(right_idx, ncol(final_df))]
  X_right = model.matrix(form, right_X, contrasts.arg = lapply(right_X[, -c(non_factor_idx, ncol(right_X))], contrasts, contrasts = FALSE))[, -1]
  in_idx_right = (1:ncol(X))[colnames(X) %in% colnames(X_right)]

  in_idx_respondent = (1:ncol(X))[-c(in_idx_left, in_idx_right)]
  #


  obs_test_stat = PO_stat(fit, in_idx_left, in_idx_right, in_idx_respondent)

  print("Initial step 1 complete: Finished computing observed test statistic. Now entering resampling step.....")

  ## parallel computing setup
  if (parallel) {
    cl <- snow::makeCluster(num_cores)
    doSNOW::registerDoSNOW(cl)
    iterations <- B
    pb <- utils::txtProgressBar(max = iterations, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    ##
    e <- foreach::foreach(j = 1:iterations, .combine = rbind,
                          .options.snow = opts) %dopar%
      {
        set.seed(j + B)
        library(CRTConjoint)
        x_df = x

        resampling_df = final_df
        a = sample(c(0,1), size =nrow(x), replace = TRUE)
        sample_idx = which(a == 1)

        kept = resampling_df[-sample_idx, ]

        b = resampling_df[sample_idx, ]
        new_y = 1 - b$Y
        new_first = b[, left_idx]
        new_second = b[, right_idx]
        name_1 = colnames(new_first)
        name_2 = colnames(new_second)
        colnames(new_first) = name_2
        colnames(new_second) = name_1
        new_df = cbind(new_second, new_first)
        new_df = cbind(new_df, b[, -c(left_idx, right_idx, ncol(b))])

        new_df$Y = new_y
        final_df_initial = rbind(kept, new_df)




        X_initial = model.matrix(form, final_df_initial, contrasts.arg = lapply(final_df_initial[, -c(non_factor_idx, ncol(final_df_initial))], contrasts, contrasts = FALSE))[, -1]

        if (!speedup) {
          best_lam = hierNet_logistic_CV(lambda, nfolds = nfolds, X = X_initial, y_var = final_df_initial$Y, tol = tol, constraint = FALSE, seed = seed)
        }

        invisible(capture.output(initial <- hierNet_logistic(as.matrix(X_initial), final_df_initial$Y, lam= best_lam, tol = tol, aa = aa)))

        resamp_TS = PO_stat(initial, in_idx_left, in_idx_right, in_idx_respondent)
        return(resamp_TS)
      }

    snow::stopCluster(cl)

  } else{
    e = vector()
    for (j in 1:B) {
      set.seed(j + B)
      source("script_AWS_source.R")
      x_df = x

      resampling_df = final_df
      a = sample(c(0,1), size =nrow(x), replace = TRUE)
      sample_idx = which(a == 1)

      kept = resampling_df[-sample_idx, ]

      b = resampling_df[sample_idx, ]
      new_y = 1 - b$Y
      new_first = b[, left_idx]
      new_second = b[, right_idx]
      name_1 = colnames(new_first)
      name_2 = colnames(new_second)
      colnames(new_first) = name_2
      colnames(new_second) = name_1
      new_df = cbind(new_second, new_first)
      new_df = cbind(new_df, b[, -c(left_idx, right_idx, ncol(b))])

      new_df$Y = new_y
      final_df_initial = rbind(kept, new_df)




      X_initial = model.matrix(form, final_df_initial, contrasts.arg = lapply(final_df_initial[, -c(non_factor_idx, ncol(final_df_initial))], contrasts, contrasts = FALSE))[, -1]

      if (!speedup) {
        best_lam = hierNet_logistic_CV(lambda, nfolds = nfolds, X = X_initial, y_var = final_df_initial$Y, tol = tol, constraint = FALSE, seed = seed)
      }
      invisible(capture.output(initial <- hierNet_logistic(as.matrix(X_initial), final_df_initial$Y, lam= best_lam, tol = tol, aa = aa)))

      e[j] = PO_stat(initial, in_idx_left, in_idx_right, in_idx_respondent)
      print(paste0("Done with task: ",j, " out of ", B, " resamples"))
    }


  }
  p_val = (length(which(e >= as.numeric(obs_test_stat[1]))) + 1)/(B + 1)

  out = list()
  out$p_val = p_val
  out$obs_test_stat = obs_test_stat
  out$resampled_test_stat = e
  out$tol = tol
  if (speedup) {
    out$lam = best_lam
  }
  out$hiernet_fit = fit
  out$seed = seed
  return(out)

}

get_carryovereffect = function(x, y, left_idx, right_idx, B, num_cores, lambda, non_factor_idx, tol, seed, parallel, profileorder_constraint, task_var, resample_func, supplyown_resamples, nfolds) {

  if (unique(sapply(x[, left_idx[!(left_idx %in% non_factor_idx)]], class)) != "factor") stop("left and right factors are not all factors please supply non_factor_idx")

  if (length(left_idx) != length(right_idx)) stop("length left idx does not match right idx")

  if((length(unique(table(x[, task_var])))) != 1) {stop("Some tasks are missing")}

  total_tasks = max(x[, task_var])
  Z_tasks = seq(2, total_tasks, by = 2)
  X_tasks = seq(1, max(Z_tasks) -1 , by = 2)

  Z_df = x[(x[, task_var] %in% Z_tasks), ][, c(left_idx, right_idx)]
  Z_left = (1:ncol(Z_df))[colnames(Z_df) %in% colnames(x)[left_idx]]
  Z_right = (1:ncol(Z_df))[!(colnames(Z_df) %in% colnames(x)[left_idx])]
  colnames(Z_df) = paste0("Z_", colnames(Z_df))
  X_df = x[(x[, task_var] %in% X_tasks), ][, c(left_idx, right_idx)]
  X_left = (1:ncol(X_df))[colnames(X_df) %in% colnames(x)[left_idx]]
  X_right = (1:ncol(X_df))[!(colnames(X_df) %in% colnames(x)[left_idx])]
  colnames(X_df) = paste0("X_", colnames(X_df))
  Y = y[(x[, task_var] %in% Z_tasks)]

  original_Xdf = X_df
  original_Zdf = Z_df

  final_df = cbind(X_df, Z_df)
  final_df$Y = Y
  just_Z = cbind(original_Zdf, Y = Y)

  if (profileorder_constraint) {
    y_new = 1 - final_df$Y
    x_df = X_df
    first = X_left
    second = X_right
    n = nrow(x_df)
    empty_df = x_df
    for (i in 1:length(first)) {
      empty_df[, first[i]] = x_df[, second[i]]
      empty_df[, second[i]] = x_df[, first[i]]
    }
    X_df = rbind(x_df, empty_df, x_df, empty_df)

    z_df = Z_df
    first = Z_left
    second = Z_right
    n = nrow(z_df)
    empty_df = z_df
    for (i in 1:length(first)) {
      empty_df[, first[i]] = z_df[, second[i]]
      empty_df[, second[i]] = z_df[, first[i]]
    }
    Z_df = rbind(Z_df, empty_df, empty_df, Z_df)

    full_Y = c(final_df$Y, y_new, y_new, final_df$Y)

    final_df = cbind(X_df, Z_df, Y = full_Y)
    just_Z = cbind(Z_df, Y = full_Y)
  }


  X = model.matrix(Y ~ ., data = final_df, contrasts.arg = lapply(final_df[, -ncol(final_df)], contrasts, contrasts = FALSE))[, -1]

  # this is because X and Z are symmetric so we can WLOG just take the first half which is always the X indexes
  idx = 1:(ncol(X)/2)

  Z_mat = model.matrix(Y ~ ., data = just_Z, contrasts.arg = lapply(just_Z[, -ncol(just_Z)], contrasts, contrasts = FALSE))[, -1]

  set.seed(seed)
  half = (nrow(Z_mat)/4)
  random_idx = suppressWarnings(split(sample(half, half, replace = FALSE), as.factor(1:nfolds)))

  if (profileorder_constraint) {
    for (i in 1:nfolds) {
      random_idx[[i]] = c(random_idx[[i]], random_idx[[i]] + half, random_idx[[i]] + 2*half, random_idx[[i]] + 3*half)
    }
  }


  best_lam = hierNet_logistic_CV(lambda, nfolds = nfolds, X = Z_mat, y_var = final_df$Y, tol = tol, constraint = profileorder_constraint, seed = seed, fold_idx = random_idx)

  print("Initial step completed: finished computing cross validated lambda")

  if (is.null(resample_func)) {
    resample_func = function(x, seed = sample(c(0, 1000), size = 1), left_idx, right_idx) {
      set.seed(seed)
      df = x[, c(left_idx, right_idx)]
      variable = colnames(x)[c(left_idx, right_idx)]
      len = length(variable)
      resampled = list()
      n = nrow(df)
      for (i in 1:len) {
        var = df[, variable[i]]
        lev = levels(var)
        resampled[[i]] = factor(sample(lev, size = n, replace = TRUE))
      }

      resampled_df = data.frame(resampled[[1]])
      for (i in 2:len) {
        resampled_df = cbind(resampled_df, resampled[[i]])
      }
      colnames(resampled_df) = colnames(df)
      return(resampled_df)
    }
  }

  chosen_initial = sample(c(0:B), size = 1)
  if (chosen_initial != 0) {
    if (is.null(supplyown_resamples)) {
      resampled_x = resample_func(x = x, left_idx = left_idx, right_idx = right_idx, seed = chosen_initial)
      resampled_x_df = resampled_x[(x[, task_var] %in% X_tasks), ]

    } else {
      resampled_x = supplyown_resamples[[chosen_initial]]
      resampled_x_df = resampled_x[(x[, task_var] %in% X_tasks), ]
    }


    colnames(resampled_x_df) = paste0("X_", colnames(resampled_x_df))

    final_df_initial = cbind(resampled_x_df, original_Zdf)
    final_df_initial$Y = Y

    if (profileorder_constraint) {
      x_df = resampled_x_df
      first = X_left
      second = X_right
      n = nrow(x_df)
      empty_df = x_df
      for (i in 1:length(first)) {
        empty_df[, first[i]] = x_df[, second[i]]
        empty_df[, second[i]] = x_df[, first[i]]
      }
      X_df = rbind(x_df, empty_df, x_df, empty_df)

      final_df_initial = cbind(X_df, Z_df, Y = full_Y)
    }


    X_initial = model.matrix(Y ~ ., data = final_df_initial, contrasts.arg = lapply(final_df_initial[, -ncol(final_df)], contrasts, contrasts = FALSE))[, -1]


    invisible(capture.output(initial <- hierNet_logistic(as.matrix(X_initial), final_df$Y, lam= best_lam, tol = tol)))

  } else {
    invisible(capture.output(initial <- hierNet_logistic(as.matrix(X), final_df$Y, lam= best_lam, tol = tol)))
  }
  aa = initial

  invisible(capture.output(fit <- hierNet_logistic(as.matrix(X), final_df$Y, lam= best_lam, tol = tol, aa = aa)))

  obs_test_stat = CO_stat(fit, idx)

  print("Initial step completed: finished computing observed test statistic. Now entering resampling step.....")

  ## parallel computing setup
  if (parallel) {
    cl <- snow::makeCluster(num_cores)
    doSNOW::registerDoSNOW(cl)
    iterations <- B
    pb <- utils::txtProgressBar(max = iterations, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    ##
    e <- foreach::foreach(j = 1:iterations, .combine = rbind,
                          .options.snow = opts) %dopar%
      {
        set.seed(j + B)
        library(CRTConjoint)

        if (is.null(supplyown_resamples)) {
          resampled_x = resample_func(x = x, left_idx = left_idx, right_idx = right_idx, seed = chosen_initial)
          resampled_x_df = resampled_x[(x[, task_var] %in% X_tasks), ]

        } else {
          resampled_x = supplyown_resamples[[j]]
          resampled_x_df = resampled_x[(x[, task_var] %in% X_tasks), ]
        }



        colnames(resampled_x_df) = paste0("X_", colnames(resampled_x_df))

        final_df_initial = cbind(resampled_x_df, original_Zdf)
        final_df_initial$Y = Y

        if (profileorder_constraint) {
          x_df = resampled_x_df
          first = X_left
          second = X_right
          n = nrow(x_df)
          empty_df = x_df
          for (i in 1:length(first)) {
            empty_df[, first[i]] = x_df[, second[i]]
            empty_df[, second[i]] = x_df[, first[i]]
          }
          X_df = rbind(x_df, empty_df, x_df, empty_df)

          final_df_initial = cbind(X_df, Z_df, Y = full_Y)
        }


        X_initial = model.matrix(Y ~ ., data = final_df_initial, contrasts.arg = lapply(final_df_initial[, -ncol(final_df_initial)], contrasts, contrasts = FALSE))[, -1]


        invisible(capture.output(initial <- hierNet_logistic(as.matrix(X_initial), final_df_initial$Y, lam= best_lam, tol = tol, aa = aa)))

        resamp_TS = CO_stat(initial, idx)

        return(resamp_TS)
      }

    snow::stopCluster(cl)

  } else {
    e = vector()
    for (j in 1:B) {

      if (is.null(supplyown_resamples)) {
        resampled_x = resample_func(x = x, left_idx = left_idx, right_idx = right_idx, seed = chosen_initial)
        resampled_x_df = resampled_x[(x[, task_var] %in% X_tasks), ]

      } else {
        resample_x = supplyown_resamples[[j]]
        resample_x_df = resampled_x[(x[, task_var] %in% X_tasks), ]
      }



      colnames(resampled_x_df) = paste0("X_", colnames(resampled_x_df))

      final_df_initial = cbind(resampled_x_df, original_Zdf)
      final_df_initial$Y = Y

      if (profileorder_constraint) {
        x_df = resampled_x_df
        first = X_left
        second = X_right
        n = nrow(x_df)
        empty_df = x_df
        for (i in 1:length(first)) {
          empty_df[, first[i]] = x_df[, second[i]]
          empty_df[, second[i]] = x_df[, first[i]]
        }
        X_df = rbind(x_df, empty_df, x_df, empty_df)

        final_df_initial = cbind(X_df, Z_df, Y = full_Y)
      }


      X_initial = model.matrix(Y ~ ., data = final_df_initial, contrasts.arg = lapply(final_df_initial[, -ncol(final_df)], contrasts, contrasts = FALSE))[, -1]


      invisible(capture.output(initial <- hierNet_logistic(as.matrix(X_initial), final_df_initial$Y, lam= best_lam, tol = tol, aa = aa)))

      e[j] = CO_stat(initial, idx)
      print(paste0("Done with task: ",j, " out of ", B, " resamples"))
    }


  }
  p_val = (length(which(e >= as.numeric(obs_test_stat[1]))) + 1)/(B + 1)

  out = list()
  out$p_val = p_val
  out$obs_test_stat = obs_test_stat
  out$resampled_test_stat = e
  out$tol = tol
  out$hiernet_fit = fit
  out$seed = seed
  out$lam = best_lam
  return(out)
}

get_fatigueeffect = function(x, y, left_idx, right_idx, B, num_cores, lambda, non_factor_idx, tol, speedup, seed, parallel, profileorder_constraint, task_var, respondent_var, nfolds) {

  if (unique(sapply(x[, left_idx[!(left_idx %in% non_factor_idx)]], class)) != "factor") stop("left and right factors are not all factors please supply non_factor_idx")

  if (length(left_idx) != length(right_idx)) stop("length left idx does not match right idx")

  if((length(unique(table(x[, task_var])))) != 1) {stop("Some tasks are missing")}

  # forcing no profile order effect constraint
  if (profileorder_constraint) {
    x_df = x
    n = nrow(x_df)
    empty_df = x_df
    y_new = 1 - y
    for (i in 1:length(left_idx)) {
      empty_df[, left_idx[i]] = x_df[, right_idx[i]]
      empty_df[, right_idx[i]] = x_df[, left_idx[i]]
    }

    final_df = rbind(x_df, empty_df)
    final_df$Y = c(y, y_new)

  } else {
    final_df = x
    final_df$Y = y
  }

  final_df = final_df[, c(left_idx, right_idx, task_var, ncol(final_df))]
  task_var_name = colnames(x)[task_var]
  final_df_var_idx = which(colnames(final_df) == task_var_name)
  set.seed(seed)

  # Creates only Z matrix to perform CV on Z
  if (speedup) {
    Z_df = final_df[, !(colnames(final_df) %in% task_var_name)]
    if (is.null(non_factor_idx)) {
      Z_long = model.matrix(Y ~ . , Z_df, contrasts.arg = lapply(Z_df[, 1:(ncol(Z_df) - 1)], contrasts, contrasts = FALSE))[, -1]
    } else {
      non_factor_names = colnames(x)[non_factor_idx]
      Z_nonfactor_idx = (1:ncol(Z_df))[(colnames(Z_df) %in% non_factor_names)]
      Z_long = model.matrix(Y ~ . , Z_df, contrasts.arg = lapply(Z_df[, -c(Z_nonfactor_idx, ncol(Z_df))], contrasts, contrasts = FALSE))[, -1]

    }

    # obtains CV lambda
    best_lam = hierNet_logistic_CV(lambda, nfolds = nfolds, X = Z_long, y_var = Z_df$Y, tol = tol, constraint = profileorder_constraint, seed = seed)

    print("Initial step completed: finished computing cross validated lambda")
  } else {
    best_lam = hierNet_logistic_CV(lambda, nfolds = nfolds, X = as.matrix(X), y_var = final_df$Y, tol = tol, constraint = profileorder_constraint, seed = seed)

  }

  form = formula(paste0("Y ~ ."))
  X = model.matrix(form, final_df, contrasts.arg = lapply(final_df[, -c(non_factor_idx, final_df_var_idx, ncol(final_df))], contrasts, contrasts = FALSE))[, -1]

  num_task = length(unique(x[, task_var]))
  respondent_idx = unique(x[, respondent_var])

  if (speedup) {
    chosen_initial = sample(c(0:B), size = 1)
    if (chosen_initial != 0) {

      # shuffle task number
      new_x = x
      for (i in 1:(length(respondent_idx))) {
        tasks = x[, task_var][(x[, respondent_var] == respondent_idx[i])]
        new_tasks = sample(tasks, replace = FALSE)
        new_x[, task_var][(x[, respondent_var] == respondent_idx[i])] = new_tasks
      }
      a = new_x[, task_var]
      resampled = final_df

      if (profileorder_constraint) {
        resampled[, (colnames(final_df) %in% task_var_name)] = c(a,a)
      } else {
        resampled[, (colnames(final_df) %in% task_var_name)] = c(a,a)
      }

      X_initial = model.matrix(form, resampled, contrasts.arg = lapply(resampled[, -c(non_factor_idx, final_df_var_idx, ncol(final_df))], contrasts, contrasts = FALSE))[, -1]

      invisible(capture.output(initial <- hierNet_logistic(as.matrix(X_initial), final_df$Y, lam= best_lam, tol = tol)))
      ## new changed way
    } else {
      invisible(capture.output(initial <- hierNet_logistic(as.matrix(X), final_df$Y, lam= best_lam, tol = tol)))
    }
    aa = initial
  } else {
    aa = NULL
  }


  invisible(capture.output(fit <- hierNet_logistic(as.matrix(X), final_df$Y, lam= best_lam, tol = tol, aa = aa)))

  I_1 = as.vector(fit$th[final_df_var_idx, ])

  I_2 = as.vector(t(fit$th[, final_df_var_idx]))

  I = (I_1 + I_2)/2

  obs_test_stat = sum(I^2)/2

  print("Initial step completed: finished computing observed test statistic. Now entering resampling step.....")

  ## parallel computing setup
  if (parallel) {
    cl <- snow::makeCluster(num_cores)
    doSNOW::registerDoSNOW(cl)
    iterations <- B
    pb <- utils::txtProgressBar(max = iterations, style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    ##
    e <- foreach::foreach(j = 1:iterations, .combine = rbind,
                          .options.snow = opts) %dopar%
      {
        set.seed(j + B)
        library(CRTConjoint)
        # shuffle task number
        new_x = x
        for (i in 1:(length(respondent_idx))) {
          tasks = x[, task_var][(x[, respondent_var] == respondent_idx[i])]
          new_tasks = sample(tasks, replace = FALSE)
          new_x[, task_var][(x[, respondent_var] == respondent_idx[i])] = new_tasks
        }
        a = new_x[, task_var]
        resampled = final_df

        if (profileorder_constraint) {
          resampled[, (colnames(final_df) %in% task_var_name)] = c(a,a)
        } else {
          resampled[, (colnames(final_df) %in% task_var_name)] = c(a,a)
        }

        X_initial = model.matrix(form, resampled, contrasts.arg = lapply(resampled[, -c(non_factor_idx, final_df_var_idx, ncol(final_df))], contrasts, contrasts = FALSE))[, -1]

        if (!speedup) {
          best_lam = hierNet_logistic_CV(lambda, nfolds = nfolds, X = as.matrix(X_initial), y_var = final_df$Y, tol = tol, constraint = profileorder_constraint, seed = seed)
        }

        invisible(capture.output(fit <- hierNet_logistic(as.matrix(X_initial), final_df$Y, lam= best_lam, tol = tol, aa = aa)))


        I_1 = as.vector(fit$th[final_df_var_idx, ])

        I_2 = as.vector(t(fit$th[, final_df_var_idx]))

        I = (I_1 + I_2)/2

        resamp_TS = sum(I^2)/2

        return(resamp_TS)
      }

    snow::stopCluster(cl)

  } else {
    e = vector()
    for (j in 1:B) {
      set.seed(j + B)
      # shuffle task number
      new_x = x
      for (i in 1:(length(respondent_idx))) {
        tasks = x[, task_var][(x[, respondent_var] == respondent_idx[i])]
        new_tasks = sample(tasks, replace = FALSE)
        new_x[, task_var][(x[, respondent_var] == respondent_idx[i])] = new_tasks
      }
      a = new_x[, task_var]
      resampled = final_df

      if (profileorder_constraint) {
        resampled[, (colnames(final_df) %in% task_var_name)] = c(a,a)
      } else {
        resampled[, (colnames(final_df) %in% task_var_name)] = c(a,a)
      }
      X_initial = model.matrix(form, resampled, contrasts.arg = lapply(resampled[, -c(non_factor_idx, final_df_var_idx, ncol(final_df))], contrasts, contrasts = FALSE))[, -1]

      if (!speedup) {
        best_lam = hierNet_logistic_CV(lambda, nfolds = nfolds, X = as.matrix(X_initial), y_var = final_df$Y, tol = tol, constraint = profileorder_constraint, seed = seed)
      }

      invisible(capture.output(fit <- hierNet_logistic(as.matrix(X_initial), final_df$Y, lam= best_lam, tol = tol, aa = aa)))


      I_1 = as.vector(fit$th[final_df_var_idx, ])

      I_2 = as.vector(t(fit$th[, final_df_var_idx]))

      I = (I_1 + I_2)/2

      e[j] = sum(I^2)/2
      print(paste0("Done with task: ",j, " out of ", B, " resamples"))
    }


  }
  p_val = (length(which(e >= as.numeric(obs_test_stat[1]))) + 1)/(B + 1)

  out = list()
  out$p_val = p_val
  out$obs_test_stat = obs_test_stat
  out$resampled_test_stat = e
  out$tol = tol
  if (speedup) {
    out$lam = best_lam
  }
  out$hiernet_fit = fit
  out$seed = seed
  return(out)
}









