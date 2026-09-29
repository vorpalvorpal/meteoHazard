# ===========================================================================
# Shared internal utilities used across the hazard / exposure functions.
# Single-sourcing the input-validation idioms and the NA-coalescing helper so
# the four hazard families validate identically and cannot silently diverge.
# ===========================================================================

# Abort (classed) if `data` is missing any of `required` columns. `arg` is the
# argument name used in the message. Returns invisibly on success.
.assert_required_cols <- function(data, required, arg = "met_data", info = NULL) {
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    msg <- c("{.arg {arg}} is missing required columns: {.val {missing_cols}}.")
    if (!is.null(info)) msg <- c(msg, "i" = info)
    cli::cli_abort(msg, class = "meteoHazard_input_error")
  }
  invisible(data)
}

# Abort (classed) if any of `cols` in `data` is non-numeric. Assumes the columns
# are present (run .assert_required_cols first). Returns invisibly on success.
.assert_numeric_cols <- function(data, cols, arg = "met_data") {
  for (col in cols) {
    if (!is.numeric(data[[col]])) {
      cli::cli_abort(
        "{.arg {arg}} column {.val {col}} must be numeric, not {.cls {class(data[[col]])}}.",
        class = "meteoHazard_input_error"
      )
    }
  }
  invisible(data)
}

# Replace NA in `x` with `default` (scalar). Vectorised; preserves length.
.na_fill <- function(x, default) {
  ifelse(is.na(x), default, x)
}


# ---- Missing-input (NA) policy --------------------------------------------- #
# Shared by dust_hazard(), litter_risk(), odour_risk(), generate_twl() and the
# layers under them: a row whose required inputs include an NA returns NA for
# that row only (the other rows are computed as normal), with ONE classed
# summary warning per call (`meteoHazard_missing_input`). Genuinely invalid,
# non-missing values (negative wind, RH > 100, ...) still error.

# Logical vector of length `n`: TRUE where any of `...` (each length 1 or `n`,
# NULL ignored) is NA at that row.
.missing_rows <- function(n, ...) {
  miss <- logical(n)
  for (x in list(...)) {
    if (is.null(x)) next
    miss <- miss | rep_len(is.na(x), n)
  }
  miss
}

# As .missing_rows(), over the `cols` of data frame `data` (absent columns are
# ignored; presence is validated by the caller).
.missing_cols_rows <- function(data, cols) {
  do.call(.missing_rows, c(list(nrow(data)), unname(lapply(cols, function(col) data[[col]]))))
}

# Emit the single summary warning for rows returned as NA because of missing
# inputs. Silent when nothing is missing. `fn` names the calling function.
.warn_missing_rows <- function(miss, fn) {
  rows <- which(miss)
  if (length(rows) == 0L) return(invisible(FALSE))
  shown <- rows[seq_len(min(10L, length(rows)))]
  more  <- if (length(rows) > length(shown)) {
    sprintf(" (and %d more)", length(rows) - length(shown))
  } else {
    ""
  }
  one <- length(rows) == 1L
  head_msg <- sprintf(
    "%d of %d row%s %s a missing input; returned NA for %s.",
    length(rows), length(miss), if (length(miss) == 1L) "" else "s",
    if (one) "has" else "have", if (one) "that row" else "those rows"
  )
  cli::cli_warn(
    c(
      paste0("{.fn {fn}}: ", head_msg),
      "i" = paste0("Row", if (one) "" else "s", ": ",
                   paste(shown, collapse = ", "), more, "."),
      "i" = "Warning class {.cls meteoHazard_missing_input} (muffle it by class to silence)."
    ),
    class = "meteoHazard_missing_input"
  )
  invisible(TRUE)
}

# Evaluate `expr` with inner missing-input warnings muffled, so a wrapper can
# emit one summary warning of its own instead of one per layer.
.muffle_missing_input <- function(expr) {
  withCallingHandlers(
    expr,
    meteoHazard_missing_input = function(w) invokeRestart("muffleWarning")
  )
}


# ---- Negative radiation components ----------------------------------------- #
# Open-Meteo's dawn / dusk hours can carry a negative diffuse_radiation beside
# a positive direct_radiation (e.g. direct 85, diffuse -62, shortwave 23 W/m^2:
# it appears to derive diffuse = shortwave - direct with direct out of step).
# Rebalance the components so the global horizontal total direct + diffuse is
# kept: a negative diffuse is set to 0 and folded into direct (floored at 0),
# then a negative direct is set to 0 and folded into diffuse (floored at 0).
# A total below -RADIATION_NEGATIVE_TOL W/m^2 is not rounding noise and still
# errors. `diffuse` may be NULL (odour without a diffuse column): a negative
# direct within the tolerance is then clamped to 0. NA is left to the
# missing-input policy. Emits one `meteoHazard_input_adjusted` message when
# anything changed. Returns list(direct, diffuse).
RADIATION_NEGATIVE_TOL <- 5

.rebalance_radiation <- function(direct, diffuse = NULL, fn,
                                 direct_name = "direct_radiation",
                                 diffuse_name = "diffuse_radiation") {
  has_dif <- !is.null(diffuse)
  dif <- if (has_dif) diffuse else rep(NA_real_, length(direct))

  total <- ifelse(is.na(direct), 0, direct) + ifelse(is.na(dif), 0, dif)
  known <- !is.na(direct) & (!has_dif | !is.na(dif))
  bad   <- known & total < -RADIATION_NEGATIVE_TOL
  if (any(bad)) {
    rows  <- which(bad)
    shown <- rows[seq_len(min(5L, length(rows)))]
    what  <- if (has_dif) {
      sprintf("row %d (%s %g + %s %g = %g)", shown, direct_name, direct[shown],
              diffuse_name, dif[shown], total[shown])
    } else {
      sprintf("row %d (%s %g)", shown, direct_name, direct[shown])
    }
    cli::cli_abort(
      c(paste0("{.fn {fn}}: radiation is clearly negative (below -",
               RADIATION_NEGATIVE_TOL, " W/m^2) in ", length(rows), " row",
               if (length(rows) == 1L) "" else "s", "."),
        "x" = paste0(paste(what, collapse = ", "),
                     if (length(rows) > length(shown)) ", ..." else "", ".")),
      class = "meteoHazard_input_error"
    )
  }

  new_dir <- direct
  new_dif <- dif
  i <- !is.na(new_dif) & new_dif < 0
  new_dir[i] <- pmax(direct[i] + new_dif[i], 0)
  new_dif[i] <- 0
  j <- !is.na(new_dir) & new_dir < 0
  new_dif[j] <- pmax(new_dif[j] + new_dir[j], 0)
  new_dir[j] <- 0
  changed <- which(i | j)

  if (length(changed) > 0L) {
    shown <- changed[seq_len(min(5L, length(changed)))]
    what  <- if (has_dif) {
      sprintf("row %d: %g + %g -> %g + %g", shown, direct[shown], dif[shown],
              new_dir[shown], new_dif[shown])
    } else {
      sprintf("row %d: %g -> %g", shown, direct[shown], new_dir[shown])
    }
    cli::cli_inform(
      c(paste0("{.fn {fn}}: ", length(changed), " hour",
               if (length(changed) == 1L) "" else "s",
               " with negative radiation adjusted",
               if (has_dif) paste0(" (", direct_name, " + ", diffuse_name,
                                   " rebalanced, total kept).")
               else paste0(" (", direct_name, " clamped to 0)."),
               ""),
        "i" = paste0(paste(what, collapse = "; "),
                     if (length(changed) > length(shown)) "; ..." else "", "."),
        "i" = "Message class {.cls meteoHazard_input_adjusted}."),
      class = "meteoHazard_input_adjusted"
    )
  }
  list(direct = new_dir, diffuse = if (has_dif) new_dif else NULL)
}
