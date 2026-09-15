# Part of nlmixr2sir.
# Importance-weight degeneracy diagnostics.
#
# Normalising on the log scale with maximum subtraction stops the weights
# overflowing, but successful normalisation says nothing about whether the
# importance sample is informative. A run in which one or two candidates carry
# almost all the weight normalises perfectly well: the retained sample is then
# a handful of repeated points dressed up as `nResample` draws, and every
# interval derived from it is far narrower than the evidence supports.
#
# These are the numbers that make that visible.

# Weights below this fraction of the uniform weight 1/n contribute effectively
# nothing. Used only for the reported count, never to alter the weights.
.sirNegligibleWeightFactor <- 0.01

# Documented warning thresholds. Deliberately loose: they are meant to catch a
# sample that is degenerate rather than merely uneven.
.sirEssFractionWarn <- 0.10
.sirMaxWeightWarn <- 0.50

#' Degeneracy summaries for a set of normalized importance weights
#'
#' @param prob Normalized resampling probabilities, one per candidate. Zero
#'   entries (failed evaluations, non-finite ratios) carry no support and are
#'   dropped rather than counted.
#' @param nSuccessful Number of candidates with a usable OFV, used as the
#'   denominator for `essFraction`.
#' @return A one-row data frame: `ess`, `essFraction`, `maxWeight`,
#'   `perplexity`, `nNonNegligible`.
#' @noRd
.sirWeightDiagnostics <- function(prob, nSuccessful = NULL) {
  p <- prob[is.finite(prob) & prob > 0]
  n <- length(p)
  if (n == 0L) {
    return(data.frame(
      ess = NA_real_,
      essFraction = NA_real_,
      maxWeight = NA_real_,
      perplexity = NA_real_,
      nNonNegligible = 0L
    ))
  }
  # Renormalise defensively: the caller's vector should already sum to one, but
  # these are diagnostics and must not depend on that.
  p <- p / sum(p)

  # Kish's effective sample size: the number of equally weighted draws that
  # would carry the same information.
  ess <- 1 / sum(p^2)
  # Perplexity, exp of the Shannon entropy: the effective number of candidates
  # actually contributing. It penalises a long tail of tiny weights less
  # harshly than ESS does, so the two together say more than either alone.
  entropy <- -sum(p * log(p))
  denom <- if (is.null(nSuccessful) || !is.finite(nSuccessful) || nSuccessful <= 0) {
    n
  } else {
    nSuccessful
  }

  data.frame(
    ess = ess,
    essFraction = ess / denom,
    maxWeight = max(p),
    perplexity = exp(entropy),
    nNonNegligible = sum(p >= .sirNegligibleWeightFactor / n)
  )
}

# Warn when the retained sample rests on too little of the proposal. Both
# thresholds are reported so the user can judge rather than just be alarmed.
.sirWarnWeightDegeneracy <- function(diag, iterNum) {
  if (!is.finite(diag$ess)) {
    return(invisible(FALSE))
  }
  essFrac <- diag$essFraction
  maxW <- diag$maxWeight
  bad_ess <- is.finite(essFrac) && essFrac < .sirEssFractionWarn
  bad_max <- is.finite(maxW) && maxW > .sirMaxWeightWarn
  if (!bad_ess && !bad_max) {
    return(invisible(FALSE))
  }

  msg <- c(
    "Iteration {iterNum}: the importance weights are concentrated on few samples."
  )
  if (bad_ess) {
    msg <- c(msg, "x" = paste0(
      "Effective sample size {round(diag$ess, 1)} is ",
      "{round(100 * essFrac, 1)}% of the usable samples ",
      "(warns below {round(100 * .sirEssFractionWarn)}%)."
    ))
  }
  if (bad_max) {
    msg <- c(msg, "x" = paste0(
      "One sample carries {round(100 * maxW, 1)}% of the weight ",
      "(warns above {round(100 * .sirMaxWeightWarn)}%)."
    ))
  }
  msg <- c(
    msg,
    "i" = "The retained sample rests on less information than its size suggests.",
    "i" = "Widen the proposal with the inflation controls, or raise {.arg nSamples}."
  )
  cli::cli_warn(msg)
  invisible(TRUE)
}
