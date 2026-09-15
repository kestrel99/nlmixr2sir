# Per-parameter diagnostics: confidence intervals by iteration, and the
# RSE/correlation matrix plot. Ported from PsN R-scripts/sir_default.R.

# Per-parameter interval for every iteration, for the proposal (all evaluated
# samples) and for the retained SIR distribution (the resampled subset). Shows
# whether the uncertainty has stopped moving between iterations.
.sirIterationIntervals <- function(x, ci = 95) {
  checkmate::assertNumber(ci, lower = 50, upper = 100)
  iterations <- attr(x, "iterations", exact = TRUE)
  if (is.null(iterations) || length(iterations) == 0L) {
    cli::cli_abort("No stored SIR iterations to summarise.")
  }
  lo <- (1 - ci / 100) / 2
  hi <- 1 - lo

  out <- lapply(seq_along(iterations), function(i) {
    raw <- iterations[[i]]$rawResults
    params <- colnames(iterations[[i]]$resampledMat)
    rbind(
      .sirIntervalFrame(
        .sirProposalRows(raw),
        params,
        i,
        "proposal",
        lo,
        hi
      ),
      .sirIntervalFrame(
        raw[raw$resamples > 0L, , drop = FALSE],
        params,
        i,
        "SIR",
        lo,
        hi
      )
    )
  })
  out <- do.call(rbind, out)
  out$asymmetry <- (out$high - out$median) / (out$median - out$low)
  rownames(out) <- NULL
  out
}

.sirIntervalFrame <- function(rows, params, iteration, type, lo, hi) {
  if (nrow(rows) == 0L) {
    return(NULL)
  }
  do.call(
    rbind,
    lapply(params, function(p) {
      v <- rows[[p]]
      v <- v[!is.na(v)]
      if (length(v) == 0L) {
        return(NULL)
      }
      q <- stats::quantile(v, probs = c(lo, 0.5, hi), names = FALSE)
      data.frame(
        iteration = as.integer(iteration),
        param = p,
        type = type,
        low = q[[1L]],
        median = q[[2L]],
        high = q[[3L]],
        stringsAsFactors = FALSE
      )
    })
  )
}

.sirIntervalPlot <- function(x, ci = 95) {
  d <- .sirIterationIntervals(x, ci = ci)
  d$group <- factor(
    paste(d$type, d$iteration),
    levels = unique(paste(d$type, d$iteration))
  )
  d$type <- factor(d$type, levels = c("proposal", "SIR"))

  ggplot2::ggplot(
    d,
    ggplot2::aes(
      x = factor(.data$iteration),
      y = .data$median,
      color = .data$type
    )
  ) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.5),
      size = 1.8
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$low, ymax = .data$high),
      position = ggplot2::position_dodge(width = 0.5),
      width = 0.3
    ) +
    ggplot2::facet_wrap(stats::as.formula("~ param"), scales = "free_y") +
    ggplot2::scale_color_manual(
      values = c(proposal = "#D94801", SIR = "#2171B5"),
      name = NULL
    ) +
    ggplot2::labs(x = "Iteration", y = paste0(ci, "% interval")) +
    ggplot2::theme_bw()
}

# The RSE/correlation matrix. Lower triangle only, RSE% on the diagonal and
# correlations off it, with the diagonal annotated by the CI asymmetry ratio
#
#   (high - median) / (median - low)
#
# which is the thing a symmetric normal-approximation covariance cannot show,
# and a large part of why SIR is run at all.
.sirRseCorData <- function(x, which = c("SIR", "proposal"), ci = 95) {
  which <- match.arg(which)
  iterations <- attr(x, "iterations", exact = TRUE)
  if (is.null(iterations) || length(iterations) == 0L) {
    cli::cli_abort("No stored SIR iterations to summarise.")
  }

  # PsN compares the first iteration's proposal against the last iteration's
  # retained SIR distribution.
  intervals <- .sirIterationIntervals(x, ci = ci)
  if (identical(which, "proposal")) {
    iter <- 1L
    raw <- iterations[[1L]]$rawResults
    rows <- .sirProposalRows(raw)
  } else {
    iter <- length(iterations)
    raw <- iterations[[iter]]$rawResults
    rows <- raw[raw$resamples > 0L, , drop = FALSE]
  }
  params <- colnames(iterations[[iter]]$resampledMat)
  mat <- as.matrix(rows[, params, drop = FALSE])
  cm <- stats::cov(mat)

  sdcor <- stats::cov2cor(cm)
  diag(sdcor) <- sqrt(diag(cm))

  estimate <- stats::setNames(x$estimate, x$param)[params]
  rse <- 100 * abs(sqrt(diag(cm)) / estimate)

  asym <- intervals[
    intervals$iteration == iter & intervals$type == which,
    ,
    drop = FALSE
  ]
  asymByParam <- stats::setNames(asym$asymmetry, asym$param)

  grid <- expand.grid(
    row = factor(params, levels = params),
    col = factor(params, levels = params),
    stringsAsFactors = FALSE
  )
  grid$value <- NA_real_
  grid$isDiagonal <- grid$row == grid$col
  grid$asymmetry <- NA_real_
  for (i in seq_len(nrow(grid))) {
    r <- as.character(grid$row[i])
    cc <- as.character(grid$col[i])
    ri <- match(r, params)
    ci_ <- match(cc, params)
    if (ri < ci_) {
      next # keep one triangle only
    }
    grid$value[i] <- if (ri == ci_) rse[[r]] else sdcor[ri, ci_]
    if (ri == ci_) {
      grid$asymmetry[i] <- asymByParam[[r]]
    }
  }
  grid <- grid[!is.na(grid$value), , drop = FALSE]
  grid$label <- round(grid$value, 2)
  # PsN's asymmetry bands
  grid$asymmetryBand <- cut(
    grid$asymmetry,
    breaks = c(-Inf, 0.5, 1, 1.25, 2, Inf),
    labels = c("<0.5", "0.5-1", "1-1.25", "1.25-2", ">2"),
    right = FALSE
  )
  grid$which <- which
  rownames(grid) <- NULL
  grid
}

.sirRseCorPlot <- function(x, which = c("SIR", "proposal"), ci = 95) {
  which <- match.arg(which)
  d <- .sirRseCorData(x, which = which, ci = ci)
  d$fill <- ifelse(d$isDiagonal, NA_real_, abs(d$value))

  ggplot2::ggplot(d, ggplot2::aes(x = .data$col, y = .data$row)) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = .data$fill),
      color = "white",
      linewidth = 0.4
    ) +
    ggplot2::geom_text(
      data = d[!d$isDiagonal, , drop = FALSE],
      ggplot2::aes(label = .data$label),
      size = 3
    ) +
    ggplot2::geom_text(
      data = d[d$isDiagonal, , drop = FALSE],
      ggplot2::aes(label = .data$label, color = .data$asymmetryBand),
      size = 3,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_gradient(
      low = "#F7FBFF",
      high = "#08519C",
      limits = c(0, 1),
      na.value = "grey92",
      name = "|correlation|"
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "<0.5" = "#A50026",
        "0.5-1" = "#F46D43",
        "1-1.25" = "#1A9850",
        "1.25-2" = "#F46D43",
        ">2" = "#A50026"
      ),
      name = "CI asymmetry",
      drop = FALSE
    ) +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = paste0(which, ": RSE (%) on the diagonal, correlation off it")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
    )
}
