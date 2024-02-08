#' Lagged Energy Distance
#'
#' Calculates pairwise lagged energy distances among multivariate time series.
#'
#' The type 1 dissimilarity between two \eqn{p}-variate time series
#' \eqn{\{X_t\}_{t=1}^{m}} and \eqn{\{Y_t\}_{t=1}^{n}} is the energy
#' distance between the lagged joint distribution upto lag \eqn{h} given by
#' \deqn{D_1(X,Y; h) = 2S_h(X,Y) - S_h(X,X) - S_h(Y,Y),}
#' where
#' \deqn{S_h(X,Y) = \sum_{s=1}^{m} \sum_{t = 1}^{n} \|X_{s:s+h} - Y_{t:t+h}\|_F^\alpha / (m-h)(n-h),}
#' \eqn{X_{t:t+h} = (X_{t},\ldots,X_{t+h})},
#' \eqn{Y_{t:t+h} = (Y_{t},\ldots,Y_{t+h})} are \eqn{p \times (h + 1)}
#' matrices, \eqn{\|\cdot\|_F} is the Frobenius norm, and \eqn{\alpha} = `a`.
#'
#' The type 2 dissimilarity is the weighted sum
#' \deqn{D_2(X,Y; w) = w_0D_1(X,Y;0)+w_1\tilde{D}(X,Y;1)+
#' w_2\tilde{D}(X,Y;2)+\ldots+w_h\tilde{D}(X,Y; h),}
#' where \eqn{w = (w_0,\ldots,w_h)} are the `weights`, and \eqn{\tilde{D}(X,Y;k)}
#' is the energy distance between the joint distributions at lag \eqn{k}.
#' The latter is given by
#' \deqn{\tilde{D}(X,Y;k) = 2T_k(X,Y) - T_k(X,X) - T_k(Y,Y),}
#' where
#' \deqn{T_k(X,Y) = \sum_{s=1}^{m} \sum_{t=1}^{n} \|X_{t.k} -Y_{t.k}\|_F^\alpha / (m-k)(n-k),}
#' and \eqn{X_{t.k} = (X_t,X_{t+k})}, \eqn{Y_{t.k} = (Y_t,Y_{t+k})} are
#' \eqn{p \times 2} matrices.
#' Note that \eqn{T_0(X,Y) = 2^{\alpha/2} S_0(X,Y)}, so that
#' \eqn{\tilde{D}(X,Y;0) = 2^{\alpha/2} D_1(X,Y;0)}.
#' @param x an object.
#' @param type character determining the type of lagged energy distance.
#' @param lag number of lags.
#' @param weights optional vector of length `lag` + 1 when `type` is 2.
#'  Overrides `lag` when incompatible.
#' @param a index for energy distance. Must be in \eqn{(0, 2]}.
#' @param ... unused.
#' @return A `dist` object.
#' @export
edist <-
    function(
        x, type = c("1", "2"), lag = 0L, weights = c(numeric(lag), 1.0), a = 1,
        ...
    ) {
        UseMethod("edist")
    }

#' @title Lagged Energy Distance for Tsibbles
#' @description See [eclust::edist].
#' @inheritParams edist
#' @param x a `tsibble` object.
#' @export
edist.tbl_ts <-
    function(x, type, lag, weights, a) {
        if (!tsibble::is_tsibble(x)) {
            rlang::abort(".data is not a tsibble object.")
        }
        mv <- tsibble::measured_vars(x)
        grouped <-
            x |>
            tsibble::group_by_key()
        sizes <-
            grouped |>
            dplyr::group_size()
        rows <-
            grouped |>
            dplyr::group_rows() |>
            unlist()
        y <-
            x[mv] |>
            as.matrix() |>
            vctrs::vec_slice(rows)
        keys <-
            grouped |>
            dplyr::group_keys()
        if (NCOL(keys) == 1) {
            keys <- keys[[1]]
        } else {
            keys <- purrr::pmap_chr(keys, paste0)
        }
        if (!is.numeric(y)) {
            rlang::abort("Measured variables are not numeric.")
        }
        lag <- as.integer(lag)
        type <- match.arg(type)
        edist_impl(y, sizes, keys, type, lag, weights, a)
    }

#' @title Lagged Energy Distance for Matrices
#' @description See [eclust::edist].
#' @inheritParams edist
#' @param x a `matrix`
#' @param sizes numeric vector of sample sizes
#' @param group_ids character vector of names for each grouped time series
#' @export
edist.matrix <-
    function(x, sizes, group_ids = NULL, type, lag, weight, a) {
        edist_impl(x, sizes, group_ids, type, lag, weights, a)
    }


edist_impl <-
    function(x, sizes, group_ids, type, lag, weights, a = 1) {
        stopifnot(
            is.numeric(x), is.matrix(x), is.numeric(sizes), is.numeric(lag)
        )
        type <- match.arg(type)
        out <-
            switch(
                type,
                "1" = energy_distance_mat(x, sizes, lag, a),
                "2" = energy_distance_mat(x, sizes, weights, a)
            )
        if (!is.null(group_ids) & (length(sizes) == length(group_ids))) {
            row.names(out) <- group_ids
        }
        return(stats::as.dist(out))
    }