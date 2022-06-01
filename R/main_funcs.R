# main functions to be exported
#' Testing whether factor matters in Conjoint Experiments
#'
#' This function takes a conjoint dataset and returns the p-value when using the
#' CRT to test if Y independent of X given Z using HierNet test statistic.
#' The function requires user to specify the outcome and all factors used in the
#' conjoint experiment and any additional respondent characteristics. By default,
#' this function assumes a uniform randomization of factor levels. In addition,
#' the function assumes the forced choice conjoint experiment and consequently assumes
#' the data to contain the left and right profile factors in separate column in the
#' dataset supplied.
#'
#' @param formula A formula object specifying the outcome variable on the left-hand side
#'  and factors of (X,Z) and respondent characteristics (V) in the right hand side.
#'  RHS variables should be separated by + signs and should only contain either left
#'  or right for each (X,Z).
#'  For example Y ~ Country_left + Education_left is sufficient as opposed to
#'  Y ~ Country_left + Country_right + Education_left + Education_right
#' @param data A dataframe containing outcome variable and all factors (X,Z,V)
#' (including both left and right profile factors). All (X,Z,V) listed in
#' the formula above are expected to be of class factor unless explicitly stated
#' in non_factor input.
#' @param X Character string specifying the variable of interest. This character
#' should match column name of X in data. For example "Country_left" is sufficient.
#' @param left Vector of column names of data that corresponds to the left profile factors
#' @param right Vector of column names of data that corresponds to the right profile factors.
#' NOTE: left and right are assumed to be the same length and the
#' order should correspond to the same variables. For example left = c("Country_left",
#' "Education_left") and right = c("Country_right", "Education_right")
#' @param design A character string of one of the following options: "Uniform",
#' "Constrained Uniform", "Nonuniform", "Manual". "Uniform" refers to a completely uniform
#' design where all (X,Z) are sampled uniformly. "Nonuniform" refers to a design where all
#' (X,Z) are sampled independently but the levels of X are not sampled uniformly.
#' If design="Nonuniform", then user should supply the non-uniform probability weights in p.
#' If in_levs is not NULL, then length of p should match in_levs. "Constrained Uniform"
#' refers to a dependent randomization design where some levels of X are only
#' possible based on certain levels of Z. See examples below. If design="Constrained Uniform"
#'  user should supply constraint_randomization list indicating the dependencies.
#' "Manual" refers to more complex conjoint designs, where the user will supply
#' their own resamples in supplyown_resamples input.
#' @param p A vector of nonuniform probability weights when design="Nonuniform".
#' Length of p should match number of levels of X or length of in_levs.
#' @param constraint_randomization List containing levels of X that can only be
#' sampled with certain values of Z.
#' The first element of constraint_randomization should contain the levels of X
#' that can only be sampled with certain values of Z, which are included in the
#' second element of the list. See example below.
#' @param supplyown_resamples List of length B that contains own resamples of X
#' when design="Manual". Each element of list should contain a dataframe
#' with the same number of rows of data and two columns for the left and right
#' profile values of X.
#' @param profileorder_constraint Boolean indicating whether to enforce profile
#' order constraint (default = TRUE)
#' @param in_levs A vector of strings that are a subset of the levels of X when
#' user wants to only test if there is any difference between a few levels of X.
#' See example below.
#' @param forced_var A character string indicating column name of Z or V that user
#' wishes to force an interaction with.
#' @param non_factor A vector of strings indicating columns of data that are not
#' factors. This should only be used for respondent characteristics (V) that are
#' not factors. For example non_factor = "Respondent_Age".
#' @param B Numeric integer value indicating the number of resamples for the CRT
#' procedure. Default value is B=200.
#' @param parallel Boolean indicating whether parallel computing should be used.
#' Default value is TRUE.
#' @param num_cores Numeric integer indicating number of cores to use when parallel=TRUE.
#' num_cores should not exceed the number of cores the user's machine can handle. Default is 2.
#' @param nfolds Numeric integer indicating number of cross-validation folds. Default is 3.
#' @param lambda Numeric vector indicating lambda used for cross-validation for
#' HierNet fit. Default lambda=c(20,30,40).
#' @param tol Numeric integer indicating acceptable tolerance for terminating optimization
#'  fit for HierNet. Default is tol=1e-3. WARNING: Do not increase as it greatly
#'  increases computation time.
#' @param speedup Boolean indicating whether to employ computational tricks to
#' make function run faster. It is always recommended to use default speedup=TRUE.
#' @param seed Seed used for CRT procedure
#' @param analysis Numeric integer indicating whether to return the top x number
#' of strongest interactions in the observed test statistic. Default analysis = 0 to not
#' return any top interactions. If analysis > 0, for example analysis = 2, then
#' the top two strongest interactions contribution to the test statistic along with which interaction is returned. NOTE: this is purely for exploratory analysis.
#'
#' @return A list containing: \item{p_val}{A numeric value for the p-value testing
#' Y independent of X given Z.}
#' \item{obs_test_stat}{A numeric value for the observed test statistic. If analysis
#'  is > 0, obs_test_stat will contain a list detailing the contribution of the main effects
#' interaction effects and the top interactions.}
#' \item{resampled_test_stat}{Matrix containing all the B resampled test statistics}
#' \item{tol}{Tolerance used for HierNet}
#' \item{lam}{Best cross-validated lambda}
#' \item{hiernet_fit}{An object of class hiernet that contains the hiernet fit for
#'  the observed test statistic}
#' \item{seed}{Seed used}
#' \item{elapsed_time}{Elapsed time}
#'
#' @export
#' @references Ham, D., Janson, L., and Imai, K. (2022)
#' Using Machine Learning to Test Causal Hypotheses in Conjoint Analysis
#' @examples
#' # Subset of Immigration Choice Conjoint Experiment Data from Hainmueller et. al. (2014).
#' data("immigrationdata")
#' form = formula("Y ~ FeatEd + FeatGender + FeatCountry + FeatReason + FeatJob +
#' FeatExp + FeatPlans + FeatTrips + FeatLang + ppage + ppeducat + ppethm + ppgender")
#' left = colnames(immigrationdata)[1:9]
#' right = colnames(immigrationdata)[10:18]
#'
#' \dontrun{
#' # Testing whether edcuation matters for immigration preferences
#' education_test = CRT_pval(formula = form, data = immigrationdata, X = "FeatEd",
#'  left = left, right = right, non_factor = "ppage", B = 100, analysis = 2)
#' education_test$p_val
#' }
#'
#' # Testing whether job matters for immigration preferences
#' constraint_randomization = list() # (Job has dependent randomization scheme)
#' constraint_randomization[["FeatJob"]] = c("Financial analyst","Computer programmer",
#' "Research scientist","Doctor")
#' constraint_randomization[["FeatEd"]] = c("Equivalent to completing two years of
#' college in the US", "Equivalent to completing a graduate degree in the US",
#'  "Equivalent to completing a college degree in the US")
#' \dontrun{
#' job_test = CRT_pval(formula = form, data = immigrationdata, X = "FeatJob",
#' left = left, right = right, design = "Constrained Uniform",
#' constraint_randomization = constraint_randomization, non_factor = "ppage", B = 100)
#' job_test$p_val
#' }
#'
#'
#' # Testing whether Mexican and European immigrants are treated indistinguishably
#' country_data = immigrationdata
#' country_data$FeatCountry = as.character(country_data$FeatCountry)
#' country_data$FeatCountry_2 = as.character(country_data$FeatCountry_2)
#' country_data$FeatCountry[country_data$FeatCountry %in% c("Germany", "France",
#' "Poland")] = "Europe"
#' country_data$FeatCountry_2[country_data$FeatCountry_2 %in% c("Germany", "France",
#'  "Poland")] = "Europe"
#' country_data$FeatCountry = factor(country_data$FeatCountry)
#' country_data$FeatCountry_2 = factor(country_data$FeatCountry_2)
#' \dontrun{
#' mexico_Europe_test = CRT_pval(formula = form, data = country_data, X = "FeatCountry",
#' left = left, right = right, design = "Nonuniform",
#' in_levs = c("Mexico", "Europe"), p = c(0.25, 0.75), non_factor = "ppage", B = 100,
#' analysis = 2)
#' }
#' \dontrun{
#' # example case with supplying own resamples
#' resample_Mexico_Europe = function(country_data) {
#'  resamples_1 = sample(c("Mexico", "Europe"), size = nrow(country_data),
#'  replace = TRUE, p = c(0.25, 0.75))
#'  resamples_2 = sample(c("Mexico", "Europe"), size = nrow(country_data),
#'  replace = TRUE, p = c(0.25, 0.75))
#'  resample_df = data.frame(resamples_1, resamples_2)
#'  return(resample_df)
#' }
#' own_resamples = list()
#' for (i in 1:100) {
#'    own_resamples[[i]] = resample_Mexico_Europe(country_data)
#' }
#' mexico_Europe_test = CRT_pval(formula = form, data = country_data, X = "FeatCountry",
#' left = left, right = right, design = "Manual",
#' in_levs = c("Mexico", "Europe"), supplyown_resamples = own_resamples,
#' non_factor = "ppage", B = 100, analysis = 2)
#' }
#' # example case with forcing with candidate gender
#' \dontrun{
#' mexico_Europe_test_force = CRT_pval(formula = form, data = country_data,
#' X = "FeatCountry", left = left, right = right, design = "Nonuniform",
#' in_levs = c("Mexico", "Europe"), p = c(0.25, 0.75), forced_var = "FeatGender",
#' non_factor = "ppage", B = 100, analysis = 0)
#' }
CRT_pval = function(formula, data, X, left, right, design = "Uniform", p = NULL, constraint_randomization = NULL, supplyown_resamples = NULL, profileorder_constraint = TRUE, in_levs = NULL, forced_var = NULL, non_factor = NULL, B = 200, parallel = TRUE,num_cores = 2,nfolds = 3, lambda = c(20, 30, 40), tol = 1e-3, speedup = TRUE, seed = sample(c(1:1000), size = 1), analysis = 0) {
  start_time = Sys.time()
  # processing stage

  if (!(design %in% c("Uniform", "Constrained Uniform", "Nonuniform", "Manual"))) stop("Design should be either Uniform, Constrained Uniform, Nonuniform, or Manual")

  left_levs = sapply(data[, left], function(x) levels(x))
  right_levs = sapply(data[, right], function(x) levels(x))

  for (i in 1:length(left_levs)) {
    check = all.equal(left_levs[[i]], right_levs[[i]])
    check = isTRUE(check)
    if (!check) stop("Left factors and right factors are not aligned")
  }

  X_right = right[left %in% X]
  X_left = left[right %in% X]

  Y_var = as.character(formula)[2]
  y = data[, Y_var]
  if (class(y) != "numeric") stop("Response is not numeric")

  X_Z_V = unlist(strsplit(as.character(formula)[3], " \\+"))
  X_Z_V = gsub(" ", "", X_Z_V)

  V = X_Z_V[!(X_Z_V %in% c(left, right))]



  x = data[, c(left, right, V)]

  xcols = (1:ncol(x))[colnames(x) %in% c(X, X_left, X_right)]
  num_x_levs = levels(data[, X])

  non_factor_idx = (1:ncol(x))[colnames(x) %in% non_factor]
  if (length(non_factor_idx) == 0) {
    non_factor_idx = NULL
  }

  left_idx = (1:ncol(x))[colnames(x) %in% left]
  right_idx = (1:ncol(x))[colnames(x) %in% right]

  if (is.null(forced_var)) {
    forced = NULL
  } else {
    if (forced_var %in% V) {
      forced = (1:ncol(x))[colnames(x) == forced_var]
    } else {
      forced_left = left[right %in% forced_var]
      forced_right = right[left %in% forced_var]
      all_forced = c(forced_var, forced_left, forced_right)
      left_forced = left_idx[colnames(x)[left_idx] %in% all_forced]
      right_forced = right_idx[colnames(x)[right_idx] %in% all_forced]
      forced = c(left_forced, right_forced)
    }
  }


  ###

  if (is.null(in_levs)) {
    resample_X = num_x_levs
  } else{
    resample_X = in_levs
  }

  if (design == "Uniform") {
    resample_func_U1 = function(resample_X, n) {
      new_samp = sample(resample_X, size = n, replace = TRUE)
      return(new_samp)
    }

    resample_func_U2 = function(resample_X, n) {
      new_samp = sample(resample_X, size = n, replace = TRUE)
      return(new_samp)
    }
    resample_func_1 = resample_func_U1; resample_func_2 = resample_func_U2
    left_allowed = right_allowed = full_X = restricted_X = NULL
  }

  if (design == "Constrained Uniform") {
    if (is.null(constraint_randomization)) stop("Please supply constraint list")

    constraint_name = names(constraint_randomization)[2]

    constraint_right = right[left %in% constraint_name]
    constraint_left = left[right %in% constraint_name]
    if (length(constraint_right) >0) {
      constraint_left = constraint_name
    } else {
      constraint_right = constraint_name
    }

    left_allowed = x[, colnames(x) %in% constraint_left] %in% constraint_randomization[[2]]
    right_allowed = x[, colnames(x) %in% constraint_right] %in% constraint_randomization[[2]]

    full_X = resample_X
    restricted_X = full_X[!(full_X %in% constraint_randomization[[1]])]

    resample_func_CU1 = function(full_X, restricted_X, n, left_allowed) {
      new_samp = sample(full_X, size = n, replace = TRUE)
      new_samp[!left_allowed] = sample(restricted_X, size = n - sum(left_allowed), replace = TRUE)
      return(new_samp)
    }

    resample_func_CU2 = function(full_X, restricted_X, n, right_allowed) {
      new_samp = sample(full_X, size = n, replace = TRUE)
      new_samp[!right_allowed] = sample(restricted_X, size = n - sum(right_allowed), replace = TRUE)
      return(new_samp)
    }
    resample_func_1 = resample_func_CU1; resample_func_2 = resample_func_CU2

  }

  if (design == "Nonuniform") {
    if (is.null(p)) stop("Please supply non-uniform probability weight p")
    resample_func_NU1 = function(resample_X, n, p, null_input = NULL) {
      new_samp = sample(resample_X, size = n, replace = TRUE, prob = p)
      return(new_samp)
    }

    resample_func_NU2 = function(resample_X, n, p, null_input = NULL) {
      new_samp = sample(resample_X, size = n, replace = TRUE, prob = p)
      return(new_samp)
    }
    resample_func_1 = resample_func_NU1; resample_func_2 = resample_func_NU2

    left_allowed = right_allowed = full_X = restricted_X = NULL

  }

  if (design == "Manual") {
    if (is.null(supplyown_resamples)) stop("Please supply own resamples")
    left_allowed = right_allowed = full_X = restricted_X = resample_func_1 = resample_func_2 = NULL
  }

  out = get_CRT_pval(x = x, y = y, xcols = xcols, left_idx = left_idx, right_idx = right_idx, design = design, B= B, forced = forced, lambda = lambda, non_factor_idx = non_factor_idx, in_levs = in_levs, analysis = analysis, p = p, resample_func_1 = resample_func_1, resample_func_2 = resample_func_2, supplyown_resamples = supplyown_resamples, resample_X = resample_X,
                     full_X = full_X, restricted_X = restricted_X, left_allowed = left_allowed, right_allowed = right_allowed,
                     tol = tol, speedup = speedup, parallel = parallel, num_cores = num_cores, profileorder_constraint = profileorder_constraint, nfolds = nfolds, seed = seed)

  end_time = Sys.time()

  out$elapsed_time = end_time - start_time

  return(out)

}


#' Testing profile order effect in Conjoint Experiments
#'
#' This function takes a conjoint dataset and returns the p-value when using the
#'  CRT to test if the profile order effect holds using HierNet test statistic.
#' The function requires user to specify the outcome and all factors used in the
#'  conjoint experiment and any additional respondent characteristics.
#' The function assumes the forced choice conjoint experiment and consequently
#' assumes the data to contain the left and right profile factors in separate column
#'  in the dataset supplied.
#'
#' @param formula A formula object specifying the outcome variable on the
#' left-hand side and factors of (X,Z) and respondent characteristics (V) in the
#'  right hand side.
#'  RHS variables should be separated by + signs and should only contain either
#'  left or right for each (X,Z).
#'  For example Y ~ Country_left + Education_left is sufficient as opposed to
#'  Y ~ Country_left + Country_right + Education_left + Education_right
#' @param data A dataframe containing outcome variable and all factors (X,Z,V)
#' (including both left and right profile factors). All (X,Z,V) listed in
#' the formula above are expected to be of class factor unless explicitly stated
#' in non_factor input.
#' @param left Vector of column names of data that corresponds to the left profile factors
#' @param right Vector of column names of data that corresponds to the right profile factors.
#'  NOTE: left and right are assumed to be the same length and the
#' order should correspond to the same variables. For example left =
#' c("Country_left", "Education_left") and right = c("Country_right", "Education_right")
#' @param non_factor A vector of strings indicating columns of data that are not factors.
#' This should only be used for respondent characteristics (V) that are not factors.
#'  For example non_factor = "Respondent_Age".
#' @param B Numeric integer value indicating the number of resamples for the CRT
#'  procedure. Default value is B=200.
#' @param parallel Boolean indicating whether parallel computing should be used.
#'  Default value is TRUE.
#' @param num_cores Numeric integer indicating number of cores to use when parallel=TRUE.
#'  num_cores should not exceed the number of cores the user's machine can handle. Default is 2.
#' @param nfolds Numeric integer indicating number of cross-validation folds. Default is 3.
#' @param lambda Numeric vector indicating lambda used for cross-validation for HierNet fit.
#' Default lambda=c(20,30,40).
#' @param tol Numeric integer indicating acceptable tolerance for terminating
#' optimization fit for HierNet. Default is tol=1e-3. WARNING: Do not increase
#' as it greatly increases computation time.
#' @param speedup Boolean indicating whether to employ computational tricks to
#' make function run faster. It is always recommended to use default speedup=TRUE.
#' @param seed Seed used for CRT procedure
#'
#' @return A list containing: \item{p_val}{A numeric value for the p-value testing
#' profile order effect.}
#' \item{obs_test_stat}{A numeric value for the observed test statistic.}
#' \item{resampled_test_stat}{Matrix containing all the B resampled test statistics}
#' \item{tol}{Tolerance used for HierNet}
#' \item{lam}{Best cross-validated lambda}
#' \item{hiernet_fit}{An object of class hiernet that contains the hiernet fit
#' for the observed test statistic}
#' \item{seed}{Seed used}
#' \item{elapsed_time}{Elapsed time}
#'
#' @export
#' @references Ham, D., Janson, L., and Imai, K. (2022)
#' Using Machine Learning to Test Causal Hypotheses in Conjoint Analysis
#' @examples
#' # Subset of Immigration Choice Conjoint Experiment Data from Hainmueller et. al. (2014).
#' data("immigrationdata")
#' form = formula("Y ~ FeatEd + FeatGender + FeatCountry + FeatReason + FeatJob +
#' FeatExp + FeatPlans + FeatTrips + FeatLang + ppage + ppeducat + ppethm + ppgender")
#' left = colnames(immigrationdata)[1:9]
#' right = colnames(immigrationdata)[10:18]
#'
#' # Testing is profile order effect is present or not in immigration data
#' \dontrun{
#' profileorder_test = CRT_profileordereffect(formula = form, data = immigrationdata,
#'  left = left, right = right, B = 100)
#' profileorder_test$p_val
#' }
CRT_profileordereffect = function(formula, data, left, right, non_factor = NULL, B = 200, parallel = TRUE,num_cores = 2,nfolds = 3, lambda = c(20, 30, 40), tol = 1e-3, speedup = TRUE, seed = sample(c(1:1000), size = 1)) {
  start_time = Sys.time()
  # processing stage

  left_levs = sapply(data[, left], function(x) levels(x))
  right_levs = sapply(data[, right], function(x) levels(x))

  for (i in 1:length(left_levs)) {
    check = all.equal(left_levs[[i]], right_levs[[i]])
    check = isTRUE(check)
    if (!check) stop("Left factors and right factors are not aligned")
  }

  Y_var = as.character(formula)[2]
  y = data[, Y_var]
  if (class(y) != "numeric") stop("Response is not numeric")

  x = data[, c(left, right)]

  non_factor_idx = (1:ncol(x))[colnames(x) %in% non_factor]
  if (length(non_factor_idx) == 0) {
    non_factor_idx = NULL
  }

  left_idx = (1:ncol(x))[colnames(x) %in% left]
  right_idx = (1:ncol(x))[colnames(x) %in% right]


  out = get_profileordereffect(x = x, y = y, left_idx = left_idx, right_idx = right_idx, B= B, lambda = lambda, non_factor_idx = non_factor_idx, tol = tol, speedup = speedup, parallel = parallel, num_cores = num_cores, nfolds = nfolds, seed = seed)

  end_time = Sys.time()

  out$elapsed_time = end_time - start_time

  return(out)

}


#' Testing carryover effect in Conjoint Experiments
#'
#' This function takes a conjoint dataset and returns the p-value when using the
#' CRT to test if the carryorder effect holds using HierNet test statistic.
#' The function requires user to specify the outcome and all factors used in the
#'  conjoint experiment and any additional respondent characteristics. By default,
#'   this function
#' assumes a uniform randomization of factor levels. The function assumes the
#' forced choice conjoint experiment and consequently assumes
#' the data to contain the left and right profile factors in separate column in
#' the dataset supplied.
#'
#' @param formula A formula object specifying the outcome variable on the left-hand
#' side and factors of (X,Z) and respondent characteristics (V) in the right hand side.
#'  RHS variables should be separated by + signs and should only contain either
#'  left or right for each (X,Z).
#'  For example Y ~ Country_left + Education_left is sufficient as opposed to
#'  Y ~ Country_left + Country_right + Education_left + Education_right
#' @param data A dataframe containing outcome variable and all factors (X,Z,V)
#' (including both left and right profile factors). All (X,Z,V) listed in
#' the formula above are expected to be of class factor unless explicitly stated
#'  in non_factor input.
#' @param left Vector of column names of data that corresponds to the left profile factors
#' @param right Vector of column names of data that corresponds to the right profile factors.
#'  NOTE: left and right are assumed to be the same length and the
#' order should correspond to the same variables. For example left =
#' c("Country_left", "Education_left") and right = c("Country_right", "Education_right")
#' @param task A character string indicating column of data that contains the task evaluation.
#'  IMPORTANT: The task variable is assumed to have no missing tasks, i.e.,
#' each respondent should have 1:J tasks. Please drop respondents with missing tasks.
#' @param design A character string of one of the following options: "Uniform" or "Manual".
#' "Uniform" refers to a completely uniform design where all (X,Z) are sampled uniformly.
#' "Manual" refers to more complex conjoint designs, where the user will supply
#' their own resamples in supplyown_resamples input.
#' @param supplyown_resamples List of length B that contains own resamples of X
#' when design="Manual". Each element of list should contain a dataframe
#' with the same number of rows of data and two columns for the left and right
#' profile values of X.
#' @param profileorder_constraint Boolean indicating whether to enforce profile
#'  order constraint (default = TRUE)
#' @param non_factor A vector of strings indicating columns of data that are not
#' factors. This should only be used for respondent characteristics (V) that are
#' not factors. For example non_factor = "Respondent_Age".
#' @param B Numeric integer value indicating the number of resamples for the CRT
#' procedure. Default value is B=200.
#' @param parallel Boolean indicating whether parallel computing should be used.
#'  Default value is TRUE.
#' @param num_cores Numeric integer indicating number of cores to use when parallel=TRUE.
#'  num_cores should not exceed the number of cores the user's machine can handle. Default is 2.
#' @param nfolds Numeric integer indicating number of cross-validation folds. Default is 3.
#' @param lambda Numeric vector indicating lambda used for cross-validation for
#' HierNet fit. Default lambda=c(20,30,40).
#' @param tol Numeric integer indicating acceptable tolerance for terminating optimization
#' fit for HierNet. Default is tol=1e-3. WARNING: Do not increase as it greatly increases
#'  computation time.
#' @param seed Seed used for CRT procedure
#'
#' @return A list containing: \item{p_val}{A numeric value for the p-value testing
#' carryover effect.}
#' \item{obs_test_stat}{A numeric value for the observed test statistic.}
#' \item{resampled_test_stat}{Matrix containing all the B resampled test statistics}
#' \item{tol}{Tolerance used for HierNet}
#' \item{lam}{Best cross-validated lambda}
#' \item{hiernet_fit}{An object of class hiernet that contains the hiernet fit for
#'  the observed test statistic}
#' \item{seed}{Seed used}
#' \item{elapsed_time}{Elapsed time}
#'
#' @export
#' @references Ham, D., Janson, L., and Imai, K. (2022)
#' Using Machine Learning to Test Causal Hypotheses in Conjoint Analysis
#' @examples
#' # Subset of Immigration Choice Conjoint Experiment Data from Hainmueller et. al. (2014).
#' data("immigrationdata")
#' form = formula("Y ~ FeatEd + FeatGender + FeatCountry + FeatReason + FeatJob +
#' FeatExp + FeatPlans + FeatTrips + FeatLang + ppage + ppeducat + ppethm + ppgender")
#' left = colnames(immigrationdata)[1:9]
#' right = colnames(immigrationdata)[10:18]
#' # Each respondent evaluated 5 tasks
#' J = 5
#' carryover_df = immigrationdata
#' carryover_df$task = rep(1:J, nrow(carryover_df)/J)
#' # Since immigration conjoint experiment had dependent randomization for several factors
#' # we supply our own resamples
#' resample_func_immigration = function(x, seed = sample(c(0, 1000), size = 1), left_idx, right_idx) {
#'  set.seed(seed)
#'  df = x[, c(left_idx, right_idx)]
#'  variable = colnames(x)[c(left_idx, right_idx)]
#'  len = length(variable)
#'  resampled = list()
#'  n = nrow(df)
#'  for (i in 1:len) {
#'    var = df[, variable[i]]
#'    lev = levels(var)
#'    resampled[[i]] = factor(sample(lev, size = n, replace = TRUE))
#'  }
#'
#'  resampled_df = data.frame(resampled[[1]])
#'  for (i in 2:len) {
#'    resampled_df = cbind(resampled_df, resampled[[i]])
#'  }
#'  colnames(resampled_df) = colnames(df)
#'
#'  #escape persecution was dependently randomized
#'  country_1 = resampled_df[, "FeatCountry"]
#'  country_2 = resampled_df[, "FeatCountry_2"]
#'  i_1 = which((country_1 == "Iraq" | country_1 == "Sudan" | country_1 == "Somalia"))
#'  i_2 = which((country_2 == "Iraq" | country_2 == "Sudan" | country_2 == "Somalia"))
#'
#'  reason_1 = resampled_df[, "FeatReason"]
#'  reason_2 = resampled_df[, "FeatReason_2"]
#'  levs = levels(reason_1)
#'  r_levs = levs[c(2,3)]
#'
#'  reason_1 = sample(r_levs, size = n, replace = TRUE)
#'
#'  reason_1[i_1] = sample(levs, size = length(i_1), replace = TRUE)
#'
#'  reason_2 = sample(r_levs, size = n, replace = TRUE)
#'
#'  reason_2[i_2] = sample(levs, size = length(i_2), replace = TRUE)
#'
#'  resampled_df[, "FeatReason"] = reason_1
#'  resampled_df[, "FeatReason_2"] = reason_2
#'
#'  #profession high skill fix
#'  educ_1 = resampled_df[, "FeatEd"]
#'  educ_2 = resampled_df[, "FeatEd_2"]
#'  i_1 = which((educ_1 == "Equivalent to completing two years of college in the US" |
#'   educ_1 == "Equivalent to completing a college degree in the US" |
#'   educ_1 == "Equivalent to completing a graduate degree in the US"))
#'  i_2 = which((educ_2 == "Equivalent to completing two years of college in the US" |
#'  educ_2 == "Equivalent to completing a college degree in the US" |
#'  educ_2 == "Equivalent to completing a graduate degree in the US"))
#'
#'
#'  job_1 = resampled_df[, "FeatJob"]
#'  job_2 = resampled_df[, "FeatJob_2"]
#'  levs = levels(job_1)
#'  # take out computer programmer, doctor, financial analyst, and research scientist
#'  r_levs = levs[-c(2,4,5, 9)]
#'
#'  job_1 = sample(r_levs, size = n, replace = TRUE)
#'
#'  job_1[i_1] = sample(levs, size = length(i_1), replace = TRUE)
#'
#'  job_2 = sample(r_levs, size = n, replace = TRUE)
#'
#'  job_2[i_2] = sample(levs, size = length(i_2), replace = TRUE)
#'
#'  resampled_df[, "FeatJob"] = job_1
#'  resampled_df[, "FeatJob_2"] = job_2
#'
#'  resampled_df[colnames(resampled_df)] = lapply(resampled_df[colnames(resampled_df)], factor )
#'
#'  return(resampled_df)
#' }
#'
#'\dontrun{
#'own_resamples = list()
#' B = 100
#' for (i in 1:B) {
#'  newdf = resample_func_immigration(carryover_df, left_idx = 1:9, right_idx = 10:18, seed = i)
#'  own_resamples[[i]] = newdf
#' }
#' carryover_test = CRT_carryovereffect(formula = form, data = carryover_df, left = left,
#' right = right, task = "task", supplyown_resamples = own_resamples, B = B)
#' carryover_test$p_val
#' }
CRT_carryovereffect = function(formula, data, left, right, task, design = "Uniform", supplyown_resamples = NULL, profileorder_constraint = TRUE, non_factor = NULL, B = 200, parallel = TRUE,num_cores = 2,nfolds = 3, lambda = c(20, 30, 40), tol = 1e-3, seed = sample(c(1:1000), size = 1)) {
  start_time = Sys.time()
  # processing stage
  if (!(design %in% c("Uniform", "Manual"))) stop("Design should be either Uniform or Manual")

  left_levs = sapply(data[, left], function(x) levels(x))
  right_levs = sapply(data[, right], function(x) levels(x))

  for (i in 1:length(left_levs)) {
    check = all.equal(left_levs[[i]], right_levs[[i]])
    check = isTRUE(check)
    if (!check) stop("Left factors and right factors are not aligned")
  }

  Y_var = as.character(formula)[2]
  y = data[, Y_var]
  if (class(y) != "numeric") stop("Response is not numeric")

  x = data[, c(left, right, task)]

  non_factor_idx = (1:ncol(x))[colnames(x) %in% non_factor]
  if (length(non_factor_idx) == 0) {
    non_factor_idx = NULL
  }

  left_idx = (1:ncol(x))[colnames(x) %in% left]
  right_idx = (1:ncol(x))[colnames(x) %in% right]

  task_var = which(colnames(x) == task)

  if (design == "Uniform") {
    resample_func = NULL
  }
  if (design == "Manual") {
    resample_func = c("Not Null")
  }
  out = get_carryovereffect(x = x, y = y, left_idx = left_idx, right_idx = right_idx, task_var = task_var, resample_func = resample_func, supplyown_resamples = supplyown_resamples, profileorder_constraint = profileorder_constraint, B= B, lambda = lambda, non_factor_idx = non_factor_idx, tol = tol, parallel = parallel, num_cores = num_cores, nfolds = nfolds, seed = seed)

  end_time = Sys.time()

  out$elapsed_time = end_time - start_time

  return(out)
}

#' Testing fatigue effect in Conjoint Experiments
#'
#' This function takes a conjoint dataset and returns the p-value when using the
#' CRT to test if the fatigue effect holds using HierNet test statistic.
#' The function requires user to specify the outcome and all factors used in the
#' conjoint experiment and any additional respondent characteristics.
#' The function assumes the forced choice conjoint experiment and consequently
#' assumes the data to contain the left and right profile factors in separate
#' column in the dataset supplied.
#'
#' @param formula A formula object specifying the outcome variable on the left-hand
#'  side and factors of (X,Z) and respondent characteristics (V) in the right hand side.
#'  RHS variables should be separated by + signs and should only contain either
#'  left or right for each (X,Z).
#'  For example Y ~ Country_left + Education_left is sufficient as opposed to
#'  Y ~ Country_left + Country_right + Education_left + Education_right
#' @param data A dataframe containing outcome variable and all factors (X,Z,V)
#' (including both left and right profile factors). All (X,Z,V) listed in
#' the formula above are expected to be of class factor unless explicitly stated
#'  in non_factor input.
#' @param left Vector of column names of data that corresponds to the left profile factors
#' @param right Vector of column names of data that corresponds to the right profile factors.
#'  NOTE: left and right are assumed to be the same length and the
#' order should correspond to the same variables. For example left =
#' c("Country_left", "Education_left") and right = c("Country_right", "Education_right")
#' @param task A character string indicating column of data that contains the task
#' evaluation. IMPORTANT: The task variable is assumed to have no missing tasks, i.e.,
#' each respondent should have 1:J tasks. Please drop respondents with missing tasks.
#' @param respondent A character string indicating column of data that contains the
#'  respondent index. The column should contain integers from 1:N indicating respondent index.
#' @param profileorder_constraint Boolean indicating whether to enforce profile order
#' constraint (default = TRUE)
#' @param non_factor A vector of strings indicating columns of data that are not factors.
#'  This should only be used for respondent characteristics (V) that are not factors.
#'   For example non_factor = "Respondent_Age".
#' @param B Numeric integer value indicating the number of resamples for the CRT procedure.
#'  Default value is B=200.
#' @param parallel Boolean indicating whether parallel computing should be used.
#' Default value is TRUE.
#' @param num_cores Numeric integer indicating number of cores to use when parallel=TRUE.
#'  num_cores should not exceed the number of cores the user's machine can handle. Default is 2.
#' @param nfolds Numeric integer indicating number of cross-validation folds. Default is 3.
#' @param lambda Numeric vector indicating lambda used for cross-validation for HierNet fit.
#' Default lambda=c(20,30,40).
#' @param tol Numeric integer indicating acceptable tolerance for terminating optimization
#' fit for HierNet. Default is tol=1e-3. WARNING: Do not increase as it greatly increases
#' computation time.
#' @param speedup Boolean indicating whether to employ computational tricks to make
#' function run faster. It is always recommended to use default speedup=TRUE.
#' @param seed Seed used for CRT procedure
#'
#' @return A list containing: \item{p_val}{A numeric value for the p-value testing
#' fatigue effect.}
#' \item{obs_test_stat}{A numeric value for the observed test statistic.}
#' \item{resampled_test_stat}{Matrix containing all the B resampled test statistics}
#' \item{tol}{Tolerance used for HierNet}
#' \item{lam}{Best cross-validated lambda}
#' \item{hiernet_fit}{An object of class hiernet that contains the hiernet fit for
#'  the observed test statistic}
#' \item{seed}{Seed used}
#' \item{elapsed_time}{Elapsed time}
#'
#' @export
#' @references Ham, D., Janson, L., and Imai, K.
#' (2022) Using Machine Learning to Test Causal Hypotheses in Conjoint Analysis
#' @examples
#' # Subset of Immigration Choice Conjoint Experiment Data from Hainmueller et. al. (2014).
#' data("immigrationdata")
#' form = formula("Y ~ FeatEd + FeatGender + FeatCountry + FeatReason + FeatJob +
#' FeatExp + FeatPlans + FeatTrips + FeatLang + ppage + ppeducat + ppethm + ppgender")
#' left = colnames(immigrationdata)[1:9]
#' right = colnames(immigrationdata)[10:18]
#' # Each respondent evaluated 5 tasks
#' J = 5
#' fatigue_df = immigrationdata
#' fatigue_df$task = rep(1:J, nrow(fatigue_df)/J)
#' fatigue_df$respondent = rep(1:(nrow(fatigue_df)/J), each = J)
#' \dontrun{
#' fatigue_test = CRT_fatigueeffect(formula = form, data = fatigue_df, left = left,
#' right = right, task = "task", respondent = "respondent", B = 100)
#' fatigue_test$p_val
#' }
CRT_fatigueeffect = function(formula, data, left, right, task, respondent, profileorder_constraint = TRUE, non_factor = NULL, B = 200, parallel = TRUE,num_cores = 2,nfolds = 3, lambda = c(20, 30, 40), tol = 1e-3, speedup = TRUE, seed = sample(c(1:1000), size = 1)) {
  start_time = Sys.time()
  # processing stage

  left_levs = sapply(data[, left], function(x) levels(x))
  right_levs = sapply(data[, right], function(x) levels(x))

  for (i in 1:length(left_levs)) {
    check = all.equal(left_levs[[i]], right_levs[[i]])
    check = isTRUE(check)
    if (!check) stop("Left factors and right factors are not aligned")
  }

  Y_var = as.character(formula)[2]
  y = data[, Y_var]
  if (class(y) != "numeric") stop("Response is not numeric")

  x = data[, c(left, right, task, respondent)]

  non_factor_idx = (1:ncol(x))[colnames(x) %in% non_factor]
  if (length(non_factor_idx) == 0) {
    non_factor_idx = NULL
  }

  left_idx = (1:ncol(x))[colnames(x) %in% left]
  right_idx = (1:ncol(x))[colnames(x) %in% right]

  task_var = which(colnames(x) == task)
  respondent_var = which(colnames(x) == respondent)

  out = get_fatigueeffect(x = x, y = y, left_idx = left_idx, right_idx = right_idx, task_var = task_var, respondent_var = respondent_var, profileorder_constraint = profileorder_constraint, B= B, lambda = lambda, non_factor_idx = non_factor_idx, tol = tol, speedup = speedup, parallel = parallel, num_cores = num_cores, nfolds = nfolds, seed = seed)

  end_time = Sys.time()

  out$elapsed_time = end_time - start_time

  return(out)

}

