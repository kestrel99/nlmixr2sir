#' nlmixr2sir package
#'
#' Sampling importance resampling tools for nlmixr2 fits.
#'
#' @keywords internal
#' @importFrom stats cov cov2cor qnorm setNames
#' @importFrom utils tail
# rxode2::ini() evaluates the OMEGA line it is handed as a lotri({...}) call in
# the *caller's* environment. lotri is an Imports of rxode2 and nlmixr2est, so
# it is never attached; without this import, OFV evaluation fails with
# "could not find function \"lotri\"" in any session where the user attached
# nlmixr2sir but not nlmixr2.
#' @importFrom lotri lotri
"_PACKAGE"
