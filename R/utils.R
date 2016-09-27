#' Repeat a data frame n number of times and return as a single data frame
#' @noRd

df_repeat <- function(x, n) {
    do.call('rbind', replicate(n, x, simplify = FALSE))
}

#' Expand one data frame to match another with no common "by" variables
#' @importFrom dplyr full_join select
#' @noRd

expand_dfs <- function(x, y) {
    x$fake___ <- 1
    y$fake___ <- 1
    combined <- full_join(x, y, by = 'fake___')
    combined <- dplyr::select(combined, -fake___)
    return(combined)
}

#' Drop a data frame column by name
#' @noRd

drop_col <- function(df, col_name) {
    out <- df[, !(names(df) %in% col_name)]
    return(out)
}
