#' Retrieve marginal means
#'
#' This function pulls emmeans::emmeans() into a pipeline that returns a dataframe
#' containing the marginal means.
#'
#' @param model model object
#' @param by variable name over which to calculate the marginal means
#'
#' @return data frame
#' @export
#'
#'

get_emmeans = function(model, by) {
  # marginal means
  model %>%
    emmeans::emmeans(specs = as.character(by)) %>%
    as.data.frame()
}
