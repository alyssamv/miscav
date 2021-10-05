#' Build a linear  model with random intercept
#'
#' This function will build a univariate LMM with random intercept
#' for a response variable that needs Box-Cox transformation. The function will
#' find the optimal lambda using the optim.bc.av function (adapted from the
#' boxcoxmix package), and produce a model using the lme4 package.
#'
#' @param f formula for linear model. Random effects may be included
#' @param data name of the dataframe being used
#'
#' @return model object
#' @export
#'
#'

lm.bc = function(f, data) {

  # search over lambda values for boxcox transform
  boxcox_lambda_search = optim.bc.av(
    y_var ~ x_var,
    groups = as.character(gp),
    data = data,
    plot.opt = 0 # only plot final
  )

  # optimal lambda for BC transform
  optimal_lambda = boxcox_lambda_search$Maximum

  ## Model with retrieved lambda
  mod = with(emmeans::make.tran("boxcox", optimal_lambda), #
             lme4::lmer(formula(f), # random intercept for subject
                        data = data,
                        na.action = na.omit,
                        REML = TRUE))

  # return the model object
  return(mod)
}
