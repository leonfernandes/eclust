#' Energy Distance
#'
#' @param .data `tsibble` object.
#' @param type character representing type of lagged energy distance.
#' @param lag number of lags.
#' @param weights optional vector of weights to use when `type` is "paired".
#' @param a coefficient for energy distance.
#' @return a `dist` object
#' @export
edist <-
    function(
        .data, type = c("upto", "paired"), lag = 0L, a = 1.0, weights = NULL
    ) {
        if (!tsibble::is_tsibble(.data)) {
            rlang::abort(".data is not a tsibble object.")
        }
        mv <- tsibble::measured_vars(.data)
        sizes <-
            .data |>
            tsibble::group_by_key() |>
            dplyr::count(name = "sizes")
        sizes <- sizes[["sizes"]]
        y <- as.matrix(.data[mv])
        if (!is.numeric(y)) {
            rlang::abort("Measured variables are not numeric.")
        }
        lag <- as.integer(lag) #always converted to integer
        type <- match.arg(type)
        out <-
            if (type == "paired") {
                if (is.null(weights)) {
                    weights <- c(numeric(lag), 1.0)
                }
                if (vctrs::vec_size(weights) != lag + 1) {
                    rlang::warn(
                        "weights and lag are incompatible. lag has been ignored."
                    )
                }
                energy_distance_mat(y, sizes, weights, a)
            } else {
                energy_distance_mat(y, sizes, lag, a)
            }
        stats::as.dist(out)
    }