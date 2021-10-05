#' Find optimal lambda for BoxCox transformation
#'
#' This function is a slightly modified version of the boxcoxmix::optim.boxcox()
#' function -- I added a progress bar. See ?boxcoxmix::optim.boxcox() for more information.
#'
#' @export
#'
#'


optim.bc.av = function (formula, groups = 1, data, K = 3, steps = 500, tol = 0.5,
                        start = "gq", EMdev.change = 1e-04, find.in.range = c(-3,
                                                                              3), s = 60, plot.opt = 3, verbose = FALSE, noformat = FALSE,
                        ...) {

  call <- match.call()
  mform <- strsplit(as.character(groups), "\\|")
  mform <- gsub(" ", "", mform)
  if (!noformat) {
    if (steps > 8)
      graphics::par(mfrow = c(4, 4), cex = 0.5)
    else graphics::par(mfrow = c(3, 3), cex = 0.5, cex.axis = 1.1)
  }
  result <- disp <- loglik <- EMconverged <- rep(0, s + 1)
  S <- 0:s
  lambda <- find.in.range[1] + (find.in.range[2] - find.in.range[1]) *
    S/s

  ## add progress bar to output
  pb <- utils::txtProgressBar(min = 0,      # Minimum value of the progress bar
                              max = (s + 1), # Maximum value of the progress bar
                              style = 3,    # Progress bar style (also available style = 1 and style = 2)
                              width = 50,   # Progress bar width. Defaults to getOption("width")
                              char = "=")   # Character used to create the bar

  for (t in 1:(s + 1)) {
    fit <- try(np.boxcoxmix(formula = formula, groups = groups,
                            data = data, K = K, lambda = lambda[t], steps = steps,
                            tol = tol, start = start, EMdev.change = EMdev.change,
                            plot.opt = plot.opt, verbose = verbose))
    if (class(fit) == "try-error") {
      cat("optim.boxcox failed using lambda=", lambda[t],
          ". Hint:  specify another range of lambda values and try again.")
      return()
    }
    EMconverged[t] <- fit$EMconverged
    result[t] <- fit$loglik
    if (!all(is.finite(result[t]))) {
      print.plot <- FALSE
    }

    # Sets the progress bar to the current state
    utils::setTxtProgressBar(pb, t)
  }
  s.max <- which.max(result)
  max.result <- result[s.max]
  lambda.max <- lambda[s.max]
  fit <- np.boxcoxmix(formula = formula, groups = groups,
                      data = data, K = K, lambda = lambda.max, steps = steps,
                      tol = tol, start = start, EMdev.change = EMdev.change,
                      plot.opt = 0, verbose = verbose)
  W <- fit$w
  P <- fit$p
  se <- fit$se
  iter <- fit$EMiteration
  names(P) <- paste("MASS", 1:K, sep = "")
  Z <- fit$mass.point
  names(Z) <- paste("MASS", 1:K, sep = "")
  Beta <- fit$beta
  Sigma <- fit$sigma
  Disp <- fit$disparity
  Disparities <- fit$Disparities
  n <- NROW(data)
  if (fit$model == "pure") {
    if (K == 1) {
      aic <- Disp + 2 * 2
      bic <- Disp + log(n) * 2
    }
    else {
      aic <- Disp + 2 * (2 * K)
      bic <- Disp + log(n) * (2 * K)
    }
  }
  else {
    if (K == 1) {
      aic <- Disp + 2 * (length(Beta) + 1)
      bic <- Disp + log(n) * (length(Beta) + 1)
    }
    else {
      aic <- Disp + 2 * (length(Beta) + 2 * K)
      bic <- Disp + log(n) * (length(Beta) + 2 * K)
    }
  }
  y <- fit$y
  yt <- fit$yt
  fitted <- fit$fitted
  fitted.transformed <- fit$fitted.transformed
  masses <- fit$masses
  ylim <- fit$ylim
  residuals <- fit$residuals
  residuals.transformed <- fit$residuals.transformed
  predicted.re <- fit$predicted.re
  Class <- fit$Class
  xx <- fit$xx
  model <- fit$model
  Disp <- fit$disparity
  Disparities <- fit$Disparities
  Loglik <- fit$loglik
  npcolors <- "green"
  ylim1 = range(result, max.result)
  maxl <- paste("Maximum profile log-likelihood:", round(max.result,
                                                         digits = 3), "at lambda=", round(lambda.max, digits = 2),
                "\n")
  if (plot.opt == 3) {
    graphics::plot(lambda, result, type = "l", xlab = expression(lambda),
                   ylab = "Profile log-likelihood", ylim = ylim1, col = "green")
    plims <- graphics::par("usr")
    y0 <- plims[3]
    graphics::segments(lambda.max, y0, lambda.max, max.result,
                       lty = 1, col = "red", lwd = 2)
    cat("Maximum profile log-likelihood:", max.result, "at lambda=",
        lambda.max, "\n")
    result <- list(call = call, y0 = y0, p = P, mass.point = Z,
                   beta = Beta, sigma = Sigma, se = se, w = W, Disparities = Disparities,
                   formula = formula, data = data, loglik = Loglik,
                   aic = aic, bic = bic, masses = masses, y = y, yt = yt,
                   All.lambda = lambda, profile.loglik = result, disparity = Disp,
                   EMconverged = EMconverged, Maximum = lambda.max,
                   mform = length(mform), ylim = ylim, fitted = fitted,
                   Class = Class, fitted.transformed = fitted.transformed,
                   predicted.re = predicted.re, residuals = residuals,
                   residuals.transformed = residuals.transformed, objective = max.result,
                   kind = 3, EMiteration = iter, ss = s, s.max = s.max,
                   npcolor = npcolors, ylim1 = ylim1, xx = xx, maxl = maxl,
                   model = model)
    class(result) <- "boxcoxmix"
  }
  else {
    result <- list(call = call, p = P, mass.point = Z, beta = Beta,
                   sigma = Sigma, se = se, w = W, Disparities = Disparities,
                   formula = formula, data = data, loglik = Loglik,
                   aic = aic, bic = bic, masses = masses, y = y, yt = yt,
                   All.lambda = lambda, profile.loglik = result, disparity = Disp,
                   EMconverged = EMconverged, Maximum = lambda.max,
                   mform = length(mform), ylim = ylim, fitted = fitted,
                   Class = Class, fitted.transformed = fitted.transformed,
                   predicted.re = predicted.re, residuals = residuals,
                   residuals.transformed = residuals.transformed, objective = max.result,
                   kind = 3, EMiteration = iter, ss = s, s.max = s.max,
                   npcolor = npcolors, ylim1 = ylim1, xx = xx, maxl = maxl,
                   model = model)
    class(result) <- "boxcoxmix"
  }
  return(result)

}
