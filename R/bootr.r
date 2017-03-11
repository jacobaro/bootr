# packages we will require
devtools::use_package("dplyr")
devtools::use_package("boot")
devtools::use_package("purrr")
devtools::use_package("cem")
devtools::use_package("parallel")
devtools::use_package("stringr")
# devtools::use_readme_md()
# devtools::document()

# imports
#' @importFrom parallel detectCores
#' @importFrom dplyr data_frame left_join %>% bind_rows
#' @importFrom boot boot boot.ci
#' @importFrom purrr by_row
#' @importFrom cem cem
#' @importFrom stringr str_detect

# get the panel of coefficients for the bootstrap
get_coefficient_panel = function(model.type, formula, data, omit, ...) {
  # run
  m = model.type(formula = formula, data = data, ...)

  # special model cases
  if(any(class(m) == "multinom")) {
    c = coef(m)
    c = paste0(rep(colnames(c), each = nrow(c)), ":", rownames(c))
  } else {
    c = names(coef(m))
  }

  # get coefficients
  if(!is.null(omit)) c = c[grep(paste(omit, collapse = "|"), c, invert = T)]

  # return the panel
  return(data_frame(coefficient = c))
}

# boot function
boot_function = function(x, i, model.type, formula, formula_match, cem_cutpoints, treatment, predictions, predictions_func, cpanel,
                         unscale, weights = 1, ...) {
  # set weights
  x$weights = weights

  # data
  dt = x[i,]

  # should we do matching?
  if(!is.null(formula_match)) {
    # weights
    dcem = cem::cem(treatment, data = as.data.frame(dt[,all.vars(formula_match)]), eval.imbalance = T, cutpoints = cem_cutpoints)
    dt$weights = dcem$w
    imb = dcem$imbalance$L1$L1
  } else {
    imb = NA
  }

  # filter the data
  dt = dplyr::filter(dt, weights > 0)

  # try to run the model
  m = NA
  try(m <- model.type(formula = formula, data = dt, weights = weights, ...), T)

  # if the model ran then go, otherwise return
  if(!is.list(m)) {
    return(m)
  }

  # get coefficients
  # special model cases
  if(any(class(m) == "multinom")) {
    vals = coef(m)
    vname = paste0(rep(colnames(vals), each = nrow(vals)), ":", rownames(vals))
    vals = as.vector(vals)
    names(vals) = vname
  } else {
    vals = coef(m)
  }

  # set coefficients
  coef = (dplyr::left_join(cpanel, dplyr::data_frame(coefficient = names(vals), value = vals), by = "coefficient"))$value
  names(coef) = cpanel$coefficient

  # should we produce predictions?
  if(is.list(predictions) & is.list(predictions$predictions)) {
    # predictions
    pr = purrr::by_row(predictions$predictions, predictions_func, model = m, full_data = dt, weighted_mean = dt$weights, .collate = "cols", .to = "value")
    pr = (pr[,stringr::str_detect(colnames(pr), "value")]) %>% t %>% as.vector
    if(!is.null(unscale)) pr = unscale(pr)
    names(pr) = paste("pr_", 1:length(pr))
  } else {
    pr = NULL
  }

  # should we produce contrasts?
  if(is.list(predictions) & is.list(predictions$contrasts) & !is.null(pr)) {
    # contrasts
    get_contrast = function(x, pr) {
      if(is.list(x)) {
        return((pr[x[[1]][1]] - pr[x[[1]][2]]) - (pr[x[[2]][1]] - pr[x[[2]][2]]))
      } else {
        return( pr[x[1]] - pr[x[2]])
      }
    }
    ct = sapply(predictions$contrasts, get_contrast, pr = pr)
    names(ct) = names(predictions$contrasts)
  } else {
    ct = NULL
  }

  # get additional model info
  s = summary(m)
  add = c("N" = length(s$residuals), "Adj-R2" = s$adj.r.squared, "AIC" = AIC(m), "Imbalance" = imb)

  # return
  return(c(coef, pr, ct, add))
}

#' Run a bootstrap
#'
#' This function is a wrapper for "boot" that runs a model and produces relevant outputs.
#' @param numruns The number of iterations.
#' @param model The statistical model (defaults to "lm").
#' @param formula The formula to use for the model.
#' @param weights Optional list of weights to use (if a matching formula is also present that will take precedence).
#' @param formula_match An optional matching formula (defaults to "NULL"). The first variable in the formual is used as the treatment.
#' @param cem_cutpoints An optional list of cutpoints for variables in the matching formula.
#' @param data The dataframe that contains all variables in 'formula' and 'formula_match'.
#' @param predictions A list with two dataframes (defaults to "NULL"): "predictions" contains a set of values to produce predictions for;
#'                    "contrasts" contains a list of vectors or lists that for which contrasts should be produced.
#' @param omit A string that contains the partial name of all variables to exclude from being returned (they will still be in the model).
#' @param unscale An optional function to unscale predictions.
#' @keywords bootstrap
#' @export
#' @examples
#' bootr()

bootr = function(numruns = 1000, model.type = lm, formula, weights = NULL, formula_match = NULL, cem_cutpoints = NULL,
                 data, ci = 0.95, predictions = NULL, omit = NULL, unscale = NULL, ...) {
  # check number of cores
  numcores = parallel::detectCores() - 1

  # generate our panel of coefficient estimates
  cpanel = get_coefficient_panel(model.type = model.type, formula = formula, data = data, omit = omit, ...)

  # get the treatment variable if we have a matching formula
  if(!is.null(formula_match)) {
    treatment = all.vars(formula_match)[1]
  } else {
    treatment = NULL
  }

  # set the weights
  if(is.null(weights)) {
    weights = 1
  }

  # run the actual bootstrap
  bout = boot(data, boot_function, numruns, parallel = "snow", ncpus = numcores,
           model.type = model.type, formula = formula, formula_match = formula_match, cem_cutpoints = cem_cutpoints, treatment = treatment,
           predictions = predictions, predictions_func = bootr.predictions, cpanel = cpanel, unscale = unscale, weights = weights, ...)

  # get results
  out = lapply(1:length(bout$t0), bootr.get_parameters, bout = bout, ci = ci) %>% bind_rows

  # identify how many runs were successfull
  out = out %>%
    bind_rows(data_frame(c_name = "Number of Replications", c_val = length(na.omit(bout$t[,1])), c_low = NA, c_high = NA, p_value = NA))

  # return the result
  return(out)
}

#' Get predictions from a model
#'
#' This function allows you to generate predictions from a model using the observed values approach.
#' @param new_data A dataframe that list the values to produce predictions at.
#' @param model The model from which predictions are produced.
#' @param full_data The full dataset used to run the model.
#' @param func The function to produce predictions (defaults to "predict")
#' @param weighted_mean A vector of observation weights (defaults to "NULL")
#' @keywords predict, bootstrap
#' @export
#' @examples
#' bootr.predictions()

bootr.predictions = function(new_data, model, full_data, func = predict, weighted_mean = NULL, .labels = T, ...) {
  # function to update data frame in line
  update.data_frame = function(fd, nd) {
    col = colSums(is.na(nd)) != nrow(nd)
    nd = nd[col == T]
    fd[, names(nd)] = nd
    return(fd)
  }

  # get predictions
  if(any(class(model) == "multinom")) {
    r = func(model, update.data_frame(full_data, new_data), type = "probs", ...)
  } else {
    r = func(model, update.data_frame(full_data, new_data), type = "response", ...)
  }

  # get the mean
  if(NCOL(r) > 1) { # for multiple outcomes
    r = colMeans(r, na.rm = T)
  } else if(is.list(r)) { # for standard errors
    if(!is.null(weighted_mean)) {
      r = c(p = weighted.mean(r[[1]], weighted_mean, na.rm = T), se = weighted.mean(r[[2]], weighted_mean, na.rm = T))
    } else {
      r = c(p = mean(r[[1]], na.rm = T), se = mean(r[[2]], na.rm = T))
    }
  } else { # otherwise
    if(!is.null(weighted_mean)) {
      r = weighted.mean(r, weighted_mean, na.rm = T)
    } else {
      r = mean(r, na.rm = T)
    }
  }

  return(r)
}

#' Get estimates from a bootstrap run.
#'
#' This function allows you to get the estimates from a bootstrap run.
#' @param i The index from the bootstrap run.
#' @param bout The return from the boostrap run.
#' @param type The type of confidence interval to generate (defaults to "perc")
#' @keywords bootstrap
#' @export
#' @examples
#' bootr.get_parameters()

bootr.get_parameters = function(i, bout, type = "perc", ci = 0.95) {
  # make sure ci is bounded
  if(ci <= 0 | ci > 1) ci = 0.95

  # get name
  name = names(bout$t0[i])
  if(is.null(name)) name = i

  # check to see if there is a problem
  if(is.na(bout$t0[i])) {
    return(data_frame("c_name" = name, "c_val" = NA, "c_low" = NA, "c_high" = NA, "p_value" = NA))
  }

  # remove infinite values
  vals = bout$t[, i]
  vals = vals[!is.infinite(vals) & !is.nan(vals) & !is.na(vals)]

  # check to make sure we have variance on the term
  if(length(unique(na.omit(vals))) < 2 | var(vals, na.rm = T) == 0) {
    return(data_frame("c_name" = name, "c_val" = bout$t0[i], "c_low" = NA, "c_high" = NA, "p_value" = NA))
  }

  # get the different confidence intervals
  p = sapply(c(0.999, 0.99, 0.95, 0.9, 0.8, ci), boot.ci, boot.out = bout, type = type, index = i)

  # make sure we have actual data
  if(is.null(p[[1]])) {
    return(NULL)
  }

  # which indices
  cl = length(p[[4, 6]]) - 1
  cu = cl + 1

  # set the p_value
  p_value = sapply(1:5, function(x) { ifelse((p[[4, x]][cl] < 0 & p[[4, x]][cu] < 0) | (p[[4, x]][cl] > 0 & p[[4, x]][cu] > 0), 1 - p[[4, x]][1], 1) }) %>% min

  # assemble our return for the 95% confidence interval
  r = data_frame("c_name" = name, "c_val" = bout$t0[i], "c_low" = p[[4, 6]][cl], "c_high" = p[[4, 6]][cu], "p_value" = p_value) # bout$t0[i] median(bout$t[,i], na.rm = T)

  # return
  return(r)
}
