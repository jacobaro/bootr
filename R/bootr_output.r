#' Write the output of a bootstrap run.
#'
#' This function allows you to write the output of a bootstrap run to an HTML file.
#' @param models The list of models to write out.
#' @keywords bootstrap
#' @export
#' @examples
#' bootr.write()

bootr.write = function(models, decimals = 3, variable_labels = NULL, additional = NULL, intercept = "(Intercept)", filename = NULL) {
  # potential formatting
  formatting = c("<HEADER>", "<LINE>")

  # get number of models
  num_models = length(models)

  # run through the models to check that they have the correct number of cols and assemble vars
  var_list = lapply(1:num_models, function(x) {
    # check column number
    if(ncol(models[[x]]) != 5) {
      stop("Model \"", names(models)[x], "\" is incorrectly formatted. Needs: variable name, coefficient, lower bound, upper bound, p-value.\n")
    }

    # return variable names
    return(models[[x]][,1])
  })

  # fix var list
  var_list = var_list %>% unlist %>% unique

  # deal with intercept
  if(!is.null(intercept)) { # because R apparently still evaluates everything in an 'if' statement
    if(intercept %in% var_list) {
      var_list = c(var_list[!(var_list %in% intercept)], intercept)
    }
  }

  # reorder variables
  vector = c(variable_labels[names(variable_labels) %in% c(var_list, formatting)], var_list[!var_list %in% names(variable_labels)])

  # correct names for those missing labels
  names(vector)[names(vector) == ""] = vector[names(vector) == ""]

  # correct null names
  if(is.null(names(vector))) names(vector) = vector

  # set variable names
  var_list = data_frame(var = names(vector), label = vector)

  # set text indent if we have a header
  if("<HEADER>" %in% names(variable_labels)) {
    p_ind = " text-indent: 15px;"
  } else {
    p_ind = ""
  }

  # set up the table header
  x = paste0("<table style='text-align:center; padding: 6px 3px; border-top: double black; border-bottom: double black'>")
  x = paste0(x, "<tr><td></td>", paste0("<td>", names(models), "<br>(", seq(1:num_models), ")</td>", collapse = ""), "</tr>")
  x = paste0(x, "<tr><td colspan='", num_models + 1, "' style='border-top: 1px solid black'></td></tr>")

  # the list of model cells
  table = list()
  add = list()

  # all the variables used
  all_vars = c()

  # go through each model
  for(i in 1:num_models) {
    # set column names
    colnames(models[[i]]) = c("c", "p", "pl", "pu", "p_value")

    # number of variables
    l = nrow(models[[i]])

    # assemble cells
    coef = models[[i]]$p %>% round(decimals) %>% format(nsmall = decimals)
    lower = models[[i]]$pl %>% round(decimals) %>% format(nsmall = decimals)
    upper = models[[i]]$pu %>% round(decimals) %>% format(nsmall = decimals)

    # set significance -- very strange problem with floating points here
    stars = case_when(models[[i]]$p_value <= 0.010001 ~ "***", models[[i]]$p_value <= 0.050001 ~ "**", models[[i]]$p_value <= 0.100001 ~ "*", T ~ "")

    # assemble coefficient and confident interval
    coef = paste0(coef, stars)
    conf_int = paste0("(", lower, ", ", upper, ")")

    # turn it into a dataframe and add column names
    cells = data_frame(models[[i]]$c, sapply(1:l, function(c) paste0("<td>", coef[c], "<br>", conf_int[c], "</td>")))
    names(cells) = c("var", i)

    # add to the full table
    table[[i]] = cells
  }

  # join into one dataframe by variable name
  for(i in 1:num_models) {
    var_list = var_list %>% left_join(table[[i]], by = "var")
  }

  # set missing to empty text
  var_list[is.na(var_list)] = "<td></td>"

  # set headers
  var_list[var_list$var %in% formatting, 3:ncol(var_list)] = ""

  # set labels
  var_list$label =
    case_when(var_list$var == "<HEADER>" ~
                paste0("<td colspan='", num_models + 1, "' style='padding-top: 6px; border-bottom: 1px solid black; text-align:left'><b>", var_list$label,"</b></td>"),
              var_list$var == "<LINE>" ~ paste0("<td colspan='", num_models + 1, "'>&nbsp;<b>", var_list$label,"</b></td>"),
              T ~ paste0("<td style='text-align:left;", p_ind, "'>", var_list$label, "</td>"))

  # collapse into a string
  var_list$var = NULL
  var_list = do.call(paste0, var_list)
  var_list = paste0("<tr>", var_list, "</tr>")
  var_list = paste0(var_list, collapse = "")

  # add data to table
  x = paste0(x, var_list)

  # add additional info
  if(!is.null(additional)) {
    # turn the additional variable into the proper format
    # comes in as a vector for each model
    # turn into a list with each variable being a vector of model values
    add = data_frame(var = lapply(additional, names) %>% dplyr::combine() %>% unique)
    for(i in 1:length(additional)) {
      add = add %>% left_join(data_frame(var = names(additional[[i]]), additional[[i]]), by = "var")
    }
    #rownames(add) = add$var
    #add$var = NULL
    add[is.na(add)] = ""

    x = paste0(x, "<tr><td colspan='", num_models + 1, "' style='border-bottom: 1px solid black'></td></tr>")

    for(i in 1:nrow(add)) {
      if(add$var[i] == "<LINE>") {
        x = paste0(x, "<tr><td colspan='", num_models + 1, "' style='border-bottom: 1px solid black'></td></tr>")
      } else {
        x = paste0(x, "<tr><td style='text-align:left'>", add$var[i], "</td>", paste0("<td>", add[i,2:ncol(add)], "</td>", collapse = ""), "</tr>")
      }

    }
  }

  # add the string and the bottom of the table
  sig_str = "Significance Level (95% CI in parentheses): * p < .1, ** p < .05, *** p < .01."
  x = paste0(x, "<tr><td colspan='", num_models + 1,
             "' style='border-top: 1px solid black; text-align:left'>",
             sig_str, "</td></tr>", "</table>")

  # write out or just return the string
  if(!is.null(filename)) {
    write(x, file = filename, sep = "")
  } else {
    return(x)
  }
}

# output bootstrap returns

#' Create summary stats for a dataframe.
#'
#' This function allows you to create summary stats for a dataframe.
#' @param df The dataframe to create summary stats for.
#' @keywords summary stats
#' @export
#' @examples
#' bootr.summary_stats()

bootr.summary_stats = function(df, variable_names) {
  # deal with factor variables
  factors = (lapply(df, class) == "factor")
  df[,factors] = lapply(df[,factors], as.integer)

  # cut NA variables
  df = df %>% select(-one_of(names(variable_names[is.na(variable_names)])))

  # create the data frame
  s = data_frame(
    "N" = sapply(df, function(x) length(x[!is.na(x)])),
    "Mean" = sapply(df, mean, na.rm = T),
    "Std. Dev." = sapply(df, sd, na.rm = T),
    "Min." = sapply(df, min, na.rm = T),
    "Max." = sapply(df, max, na.rm = T))

  # round
  s = s %>% round(3)

  # set variable name
  s$var = colnames(df)

  # reorder and name

  # prune variables not actually present
  variable_names = variable_names[names(variable_names) %in% s$var]

  # get remaining variables not named
  remaining = s$var[!(s$var %in% names(variable_names))]

  # reorder to match variable names
  s = s %>% slice(match(c(names(variable_names), remaining), var))

  # set the name
  s = s %>% left_join(data_frame(Name = variable_names, var = names(variable_names)), by = "var")
  s$Name[is.na(s$Name)] = s$var[is.na(s$Name)]

  # reorder columns
  s = s %>% select(Name, everything()) %>% select(-var)

  # return
  return(s)
}

# turn an lm into a format that is printable here
lm_to_bs = function(m, summary = T) {
  if(summary) {
    # get summary
    sm = summary(m)
  } else {
    sm = m[1:nrow(m),1:4] %>% as.data.frame
    colnames(sm) = c("p.coeff", "se", "t", "p.pv")
    rownames(sm) = rownames(m)
  }

  # create df
  df = data_frame(c = rownames(sm), p = sm$p.coeff, pl = sm$p.coeff-qnorm(0.975)*sm$se, pu = sm$p.coeff+qnorm(0.975)*sm$se, p_value = sm$p.pv)

  # return
  return(df)
}

# pull out additional reporting data for a mixed effects model
bootr.out.add.lme = function(m, names) {
  # get the RE variance
  t = as.data.frame(VarCorr(m))

  # set the names
  t$grp[!is.na(match(t$grp, names(names)))] = names[!is.na(match(names(names), t$grp))]
  t$var1[!is.na(match(t$var1, names(names)))] = names[!is.na(match(names(names), t$var1))]
  t$var21[!is.na(match(t$var2, names(names)))] = names[!is.na(match(names(names), t$var2))]

  # combine
  t = t %>% unite(var_name, grp, var1, var2, sep = ":")
  t$var_name = t$var_name %>% str_replace_all(":NA", "")

  # set raneffects
  ranef = paste0(round(t$vcov, 3), " (", round(t$sdcor, 3), ")")
  names(ranef) = t$var_name

  # return
  return(c("Random Effects:" = NA, ranef, "<LINE>" = NA, "N" = nrow(m@frame), "AIC" = round(AIC(m), 3)))
}

# pull out additional reporting data for a linear model
bootr.out.add.lm = function(m) {
  # get summary
  ms = summary(m)

  # get r squared
  if(!is.null(ms$r.sq)) {
    ars = ms$r.sq # for GAM
  } else if(!is.null(ms$adj.r.squared)) {
    ars = ms$adj.r.squared # for LM
  } else {
    ars = 1-(m$deviance/m$null.deviance) # assemble adjusted r2 by comparing l1 to l0
  }

  # return
  return(c("N" = length(m$residuals), "AIC" = round(AIC(m), 3), "Adj R2" = round(ars, 3)))
}
