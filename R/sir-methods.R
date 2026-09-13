# Part of nlmixr2sir. Split out of the original single-file R/sir.R.
# S3 print and plot methods for nlmixr2SIR.

.sirParameterPlotData <- function(x) {
  resampled <- attr(x, "resampledMat")
  if (is.null(resampled)) {
    cli::cli_abort("SIR object is missing {.field resampledMat}.")
  }
  data.frame(
    param = rep(colnames(resampled), each = nrow(resampled)),
    value = as.vector(resampled),
    stringsAsFactors = FALSE
  )
}

.sirRawResults <- function(x) {
  iterations <- attr(x, "iterations")
  if (is.null(iterations) || length(iterations) == 0L) {
    cli::cli_abort("SIR object is missing iteration raw results.")
  }
  out <- lapply(seq_along(iterations), function(i) {
    raw <- iterations[[i]]$rawResults
    raw$iter <- iterations[[i]]$iterSummary$iter
    raw
  })
  do.call(rbind, out)
}

.sirSignifDataFrame <- function(x, digits) {
  out <- as.data.frame(x)
  num_cols <- vapply(out, is.numeric, logical(1L))
  out[num_cols] <- lapply(out[num_cols], signif, digits = digits)
  out
}

#' Print an SIR result
#'
#' Prints the final parameter uncertainty summary and, when available, the
#' per-iteration SIR diagnostics stored on the `nlmixr2SIR` object returned by
#' [runSIR()].
#'
#' @param x An object returned by [runSIR()].
#' @param ... Unused.
#' @param digits Number of significant digits to print.
#' @return Invisibly returns `x`.
#' @export
print.nlmixr2SIR <- function(x, ..., digits = 3) {
  checkmate::assertDataFrame(x)
  checkmate::assertCount(digits, positive = TRUE)

  output_dir <- attr(x, "outputDir", exact = TRUE)
  cli::cli_h1("SIR Summary")
  if (!is.null(output_dir)) {
    cli::cli_alert_info("Output directory: {.path {output_dir}}")
  }

  summary_cols <- intersect(
    c("param", "estimate", "sd", "rse", "p2.5", "p50", "p97.5"),
    names(x)
  )
  cli::cli_h2("Final parameter uncertainty")
  print(
    .sirSignifDataFrame(x[, summary_cols, drop = FALSE], digits = digits),
    row.names = FALSE
  )

  iter_summary <- attr(x, "iterationSummary", exact = TRUE)
  if (!is.null(iter_summary) && nrow(iter_summary) > 0L) {
    iter_cols <- intersect(
      c(
        "iter",
        "nSamples",
        "nAttempted",
        "nSuccessful",
        "nResample",
        "nResampled",
        "nNegativeDOFV",
        "minDOFV"
      ),
      names(iter_summary)
    )
    cli::cli_h2("Iteration diagnostics")
    print(
      .sirSignifDataFrame(
        iter_summary[, iter_cols, drop = FALSE],
        digits = digits
      ),
      row.names = FALSE
    )
  }

  invisible(x)
}

#' Plot an SIR result
#'
#' Plots diagnostics for an `nlmixr2SIR` object. The default plot shows the
#' final resampled parameter distributions with reference estimates marked.
#' Additional diagnostic views show dOFV distributions or resampling
#' probabilities by iteration.
#'
#' @param x An object returned by [runSIR()].
#' @param y Unused; included for S3 compatibility.
#' @param type Plot type: `"parameters"`, `"dofv"`, or `"resampling"`.
#' @param bins Number of histogram bins for parameter and dOFV plots.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
#' @importFrom ggplot2 .data
plot.nlmixr2SIR <- function(
  x,
  y,
  ...,
  type = c("parameters", "dofv", "resampling"),
  bins = 30
) {
  type <- match.arg(type)
  checkmate::assertCount(bins, positive = TRUE)

  if (type == "parameters") {
    plot_df <- .sirParameterPlotData(x)
    estimate_df <- data.frame(
      param = x$param,
      estimate = x$estimate,
      stringsAsFactors = FALSE
    )
    return(
      ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$value)) +
        ggplot2::geom_histogram(
          bins = bins,
          fill = "#6BAED6",
          color = "white"
        ) +
        ggplot2::geom_vline(
          data = estimate_df,
          ggplot2::aes(xintercept = .data$estimate),
          color = "#D94801",
          linewidth = 0.6
        ) +
        ggplot2::facet_wrap(stats::as.formula("~ param"), scales = "free") +
        ggplot2::labs(x = "Parameter value", y = "Resampled vectors") +
        ggplot2::theme_bw()
    )
  }

  raw_df <- .sirRawResults(x)
  raw_df <- raw_df[raw_df$sample_id > 0, , drop = FALSE]

  if (type == "dofv") {
    return(
      ggplot2::ggplot(raw_df, ggplot2::aes(x = .data$dOFV)) +
        ggplot2::geom_histogram(
          bins = bins,
          fill = "#74C476",
          color = "white"
        ) +
        ggplot2::facet_wrap(stats::as.formula("~ iter"), scales = "free_y") +
        ggplot2::labs(x = "dOFV", y = "Sampled vectors") +
        ggplot2::theme_bw()
    )
  }

  ggplot2::ggplot(
    raw_df,
    ggplot2::aes(x = .data$sample_id, y = .data$probability_resample)
  ) +
    ggplot2::geom_col(fill = "#9E9AC8") +
    ggplot2::facet_wrap(stats::as.formula("~ iter"), scales = "free_x") +
    ggplot2::labs(x = "Sample id", y = "Probability resample") +
    ggplot2::theme_bw()
}
