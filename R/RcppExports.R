# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

energy_distance_xy <- function(x, y, h = 0L, a = 1.0) {
    .Call(`_eclust_energy_distance_xy`, x, y, h, a)
}

energy_distance_mat <- function(mat, sizes, h, a = 1.0) {
    .Call(`_eclust_energy_distance_mat`, mat, sizes, h, a)
}

