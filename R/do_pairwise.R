#' Do pairwise comparisons from a model object
#'
#' This function pulls emmeans::emmeans() and emmeans::pairs() to do pairwise
#' comparisons with user-chosen multiple testing correction method. In addition
#' to providing the means and p-values, it also produces confidence intervals based
#' on the standard error.
#'
#' @param model model object
#' @param by variable name over which to calculate the marginal means
#' @param correction correction method. (e.g. "bonferroni", "holm", "BH")
#' @param ci width of confidence interval. Defaults to 0.95.
#'
#' @return data frame of pairwise comparisons, adjusted p-values, and confidence intervals.
#' @export
#'
#'


do_pairwise = function(model, by, correction, ci = 0.95) {
  z = abs(qnorm( (1 - ci)/2 )) # z-value for CI calculation

  model %>%
    emmeans(specs = as.character(by)) %>%
    pairs(adjust = as.character(correction)) %>%
    as.data.frame() %>%
    mutate(p.value = round(p.value, 4),
           lower.CL = round(estimate - z*SE, 2),
           upper.CL = round(estimate + z*SE, 2)) %>%
    dplyr::select(-df, -t.ratio)

  }
