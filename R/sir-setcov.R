# Registering the SIR covariance with the fit, so that
# nlmixr2est::setCov(fit, "sir") switches the fit's reported uncertainty to the
# SIR result. Mirrors nlmixr2boot's .registerBootCovList() /
# .bootstrapCovAsFitCov().
#
# nlmixr2boot has to guess how its parameter names correspond to fit$cov
# rownames, reconstructing "cov.<eta1>.<eta2>" and "om.<eta>" by pattern. That
# guesswork is unnecessary here: .sirParamSpace() already carries the exact
# sirName <-> covName correspondence, so the mapping is a lookup.

# Reshape the SIR covariance into fit$cov's names and order, or return NULL
# with an explanation if it cannot be done safely.
.sirCovAsFitCov <- function(fit, covSir, ps = .sirParamSpace(fit)) {
  fitCov <- fit$cov
  if (!is.matrix(fitCov) || is.null(rownames(fitCov))) {
    cli::cli_inform(c(
      "i" = "No {.code fit$cov} to match; skipping {.fn setCov} registration."
    ))
    return(NULL)
  }
  if (!is.matrix(covSir) || is.null(rownames(covSir))) {
    return(NULL)
  }

  fitNames <- rownames(fitCov)
  mapped <- ps$covName[match(rownames(covSir), ps$sirName)]

  if (anyNA(mapped) || !setequal(mapped, fitNames)) {
    cli::cli_inform(c(
      "i" = "SIR parameters do not match {.code fit$cov} exactly; skipping {.fn setCov} registration.",
      "i" = "This is expected when SIR ran on a subset, or on a fit whose covariance step failed."
    ))
    return(NULL)
  }

  out <- covSir
  dimnames(out) <- list(mapped, mapped)
  out <- out[fitNames, fitNames, drop = FALSE]

  if (inherits(try(chol(out), silent = TRUE), "try-error")) {
    cli::cli_inform(c(
      "i" = "The SIR covariance is not positive definite; skipping {.fn setCov} registration."
    ))
    return(NULL)
  }
  out
}

# Merge a covariance into fit$env$covList under `label`, leaving any other
# registered matrices in place.
.sirRegisterCovList <- function(fit, label, covMat) {
  if (is.null(covMat)) {
    return(invisible(FALSE))
  }
  env <- tryCatch(fit$env, error = function(e) NULL)
  if (!is.environment(env)) {
    return(invisible(FALSE))
  }
  existing <- if (exists("covList", envir = env, inherits = FALSE)) {
    get("covList", envir = env)
  } else {
    list()
  }
  existing[[label]] <- covMat
  assign("covList", existing, envir = env)
  invisible(TRUE)
}

# Called at the end of a run: register the empirical SIR covariance as "sir".
.sirRegisterCov <- function(
  fit,
  summary,
  ps = .sirParamSpace(fit),
  label = "sir"
) {
  covMat <- attr(summary, "covMatrix", exact = TRUE)
  if (is.null(covMat)) {
    return(invisible(FALSE))
  }
  .sirRegisterCovList(fit, label, .sirCovAsFitCov(fit, covMat, ps))
}
