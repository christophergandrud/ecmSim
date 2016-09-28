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

#' Convert \code{ci} interval from percent to proportion and check if valid
#' @noRd

ci_check <- function(x) {
    if (x > 1 & x <= 100) x <- x / 100
    if (x <= 0 | x > 1) {
        stop(sprintf("%s is not a valid central interval.", x),
             call. = FALSE)
    }
    return(x)
}

#' Constrict a data frame of simulated values to a central interval
#' @param sims_scenarios a data frame of simulated quantities of interest and
#' a column grouping them by fitted scenario.
#' @param scenario_var character string of the variable name marking the
#' scenarios.
#' @param qi_var character string of the name of the variable with the
#' simulated quantity of interest values.
#' @param ci numeric value indicating the central interval. Must be in (0, 1].
#'
#' @importFrom dplyr bind_rows
#' @noRd

qi_central_interval <- function(sims_scenarios, scenario_var = 'scenario_',
                                qi_var = 'qi_', ci = 0.95)
{
    qi_ <- NULL
    lower <- (1 - ci)/2
    upper <- 1 - lower

    names(sims_scenarios)[names(sims_scenarios) == qi_var] <- 'qi_'

    qi_list <- split(sims_scenarios, sims_scenarios[[scenario_var]])
    qi_list <- lapply(seq_along(qi_list), function(x){
        lower_bound <- quantile(qi_list[[x]][, 'qi_'], prob = lower)
        upper_bound <- quantile(qi_list[[x]][, 'qi_'], prob = upper)
        subset(qi_list[[x]], qi_ >= lower_bound & qi_ <= upper_bound)
    })

    out <- data.frame(bind_rows(qi_list))
    names(out)[names(out) == 'qi_'] <- qi_var

    return(out)
}
