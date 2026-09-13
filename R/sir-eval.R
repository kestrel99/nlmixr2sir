# Part of nlmixr2sir. Split out of the original single-file R/sir.R.
# OFV evaluation of sampled parameter vectors.

# Step 4 -----------------------------------------------------------------------

#' Evaluate OFV for a matrix of sampled parameter vectors
#'
#' For each row of `paramSamples`, seeds the model's ini block with the
#' supplied values (theta and/or lower-triangle omega elements) and calls
#' `nlmixr2(est = "focei")` with `maxOuterIterations = 0`, making it the
#' nlmixr2 equivalent of NONMEM `MAXEVAL=0`.  THETA column names must match
#' parameter names in `fit$iniDf`; omega columns use the SIR lower-triangle
#' proposal names.
#'
#' Parameters NOT present as columns in `paramSamples` retain their estimated
#' values from `fit`.
#'
#' @param fit An nlmixr2 fit object (carries the data internally).
#' @param paramSamples Numeric matrix, one row per sample.  Column names may
#'   include THETA names from `fit$iniDf` and SIR omega lower-triangle names.
#' @param workers Passed to `.withWorkerPlan()`: `NULL` (keep current plan),
#'   `1` (force sequential), a positive integer, or `"auto"`.
#' @param rxThreads Integer, `"auto"`, or `NULL`; rxode2 OpenMP threads per
#'   worker. Required by `.withWorkerPlan()` whenever `workers > 1`, since
#'   each worker runs its own thread pool.
#' @return Named numeric vector of length `nrow(paramSamples)`.  Entries are
#'   `NA_real_` for rows that produced an error during evaluation.
#' @noRd
sirEvalOFV <- function(fit, paramSamples, workers = NULL, rxThreads = NULL) {
  checkmate::assertClass(fit, "nlmixr2FitCore")
  checkmate::assertMatrix(
    paramSamples,
    mode = "numeric",
    min.rows = 1L,
    min.cols = 1L
  )

  ps <- .sirParamSpace(fit)
  col_names <- colnames(paramSamples)
  theta_cols <- intersect(
    col_names,
    ps$sirName[ps$kind %in% c("theta", "sigma")]
  )
  omega_cols <- intersect(
    col_names,
    ps$sirName[ps$kind %in% c("omegaDiag", "omegaOffdiag")]
  )
  has_omega <- length(omega_cols) > 0L

  if (length(theta_cols) == 0L && !has_omega) {
    cli::cli_abort(
      "No column names in {.arg paramSamples} match any parameter in {.arg fit}."
    )
  }

  base_omega <- fit$omega
  base_theta <- fit$theta

  eval_one <- function(i) {
    row <- paramSamples[i, ]

    theta_vals <- base_theta
    if (length(theta_cols) > 0L) {
      theta_vals[theta_cols] <- row[theta_cols]
    }

    ini_args <- as.list(theta_vals)

    if (has_omega) {
      omega_mat <- .sirReconstructOmega(ps, row, base_omega)
      eta_names <- rownames(omega_mat)
      lt_vals <- unlist(
        lapply(seq_len(nrow(omega_mat)), function(r) {
          omega_mat[r, seq_len(r)]
        }),
        use.names = FALSE
      )
      lhs <- paste(eta_names, collapse = " + ")
      lt_txt <- format(lt_vals, scientific = TRUE, digits = 17, trim = TRUE)
      rhs <- if (length(lt_vals) == 1L) {
        lt_txt
      } else {
        paste0("c(", paste(lt_txt, collapse = ", "), ")")
      }
      omega_expr <- str2lang(paste(lhs, "~", rhs))
      ini_args <- c(ini_args, list(omega_expr))
    }

    tryCatch(
      {
        model_new <- suppressMessages(
          do.call(rxode2::ini, c(list(x = fit), ini_args))
        )
        f <- suppressMessages(
          nlmixr2est::nlmixr2(
            model_new,
            est = "focei",
            control = nlmixr2est::foceiControl(
              calcTables = FALSE,
              covMethod = "",
              compress = FALSE,
              maxOuterIterations = 0L,
              print = 0L
            )
          )
        )
        list(objf = f$objf, error = NA_character_)
      },
      # The message is kept rather than discarded: a configuration problem
      # (a missing import, an unloadable model) fails every sample
      # identically, and reporting only "all evaluations failed" hides why.
      error = function(e) {
        list(objf = NA_real_, error = conditionMessage(e))
      }
    )
  }

  results <- nlmixr2utils::.withWorkerPlan(
    workers,
    rxThreads = nlmixr2utils::resolveRxThreads(workers, rxThreads),
    {
      # nolint: object_usage_linter.
      nlmixr2utils::.plap(
        # nolint: object_usage_linter.
        seq_len(nrow(paramSamples)),
        eval_one,
        .label = function(i) sprintf("sample %d", i)
      )
    }
  )

  objf <- vapply(results, function(r) r$objf, numeric(1L))
  errs <- unique(stats::na.omit(vapply(
    results,
    function(r) r$error,
    character(1L)
  )))
  if (length(errs) > 0L) {
    attr(objf, "evalErrors") <- as.character(errs)
  }
  objf
}
