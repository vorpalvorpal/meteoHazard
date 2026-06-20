#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr case_when
#' @importFrom purrr pmap_dbl
#' @importFrom cli cli_h1 cli_h2 cli_alert cli_alert_info cli_alert_success
#'   cli_alert_warning cli_alert_danger cli_ul cli_progress_bar
#'   cli_progress_update cli_progress_done
#' @import S7
## usethis namespace: end
NULL

.onLoad <- function(libname, pkgname) {
  S7::methods_register()
}
