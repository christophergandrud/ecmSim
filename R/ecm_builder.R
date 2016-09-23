#' Create simulations for long-term effects from error correction models (ECM)
#'
#' @param obj a fitted model object for an ECM model
#' @param baseline_df a data frame with fitted values for the baseline scenario
#' before the shock to \code{lag_dv}. Note column names should match coefficient
#' names in \code{obj}. Variables not include are assumed to be 0.
#' Note also that change variables should be 0.
#' @param lag_dv charachter string identifying the name of the lagged dependent
#' variable in \code{obj}.
#' @param lag_iv charachter string identifying the name of the lagged
#' independent variable in \code{obj} that experiences a shock.
#' @param d_iv charachter string identifying the name of the change in
#' \code{lag_dv} between the current and previous time, i.e. the size of the
#' one period shock to \code{lag_dv}.
#' @param iv_shock numeric shock to \code{iv}.
#' @param t_extent numeric specifying the time points from the shock to
#' simulate the long-term effects of the shock for. The default is 5 time
#' points.
#' @param nsims numeric. Number of simulations to draw.
#' @param mu an optional vector giving the means of the variables estimated
#' from an ECM model. If \code{obj} is supplied then \code{mu} is ignored.
#' @param Sigma an optional positive-definite symmetric matrix estimated
#' from an ECM model. The matrix specifies the covariance matrix of the
#' variables. If \code{obj}  is supplied then \code{Sigma} is ignored.
#'
#' @importFrom coreSim b_sim qi_builder
#' @export

ecm_builder <- function(obj, baseline_df, lag_dv,
                        lag_iv, d_iv, iv_shock, t_extent = 5,
                        nsim = 1000,
                        mu, Sigma)
{
    if (t_extent < 3) {
        message('t_extent must be 3 or more time points.\nForcing t_extent = 3.')
        t_extent = 3
    }

    if (missing(lag_dv)) {
        message('lag_dv not supplied. Assuming first column of baseline_df is the lagged dependent variable.\n')
        lag_dv <- names(baseline_df)[1]
    }

    # Baseline cenario
    baseline_scenario <- df_repeat(baseline_df, n = t_extent)
    baseline_scenario$time__ <- 1:t_extent
    return(baseline_scenario)
    baseline_scenario[, lag_dv][baseline_scenario$time__ > 1] <- NA

    # Create shock fitted values
    shocked <- df_repeat(baseline_df, n = t_extent)
    shocked$time__ <- 1:t_extent
    shocked[2, d_iv] <- iv_shock
    shocked[is.na(shocked)] <- 0
    shocked[3:t_extent, lag_iv] <- shocked[1, lag_iv] + iv_shock
    shocked[2:t_extent, lag_dv] <- NA

    scenarios <- list(baseline = baseline_scenario,
                      shocked = shocked)

    # Simulate parameters
    if (!missing(obj))
        param_sims <- b_sim(obj = obj, nsim = nsim)
    else if (!missing(mu) & !missing(Sigma))
        param_sims <- b_sim(mu = mu, Sigma = Sigma, nsim = nsim)


    sims <- data.frame()
    sims <- lapply(seq_along(scenarios), function(x) {
        temp_scen <- scenarios[[x]]

        temp_sims <- data.frame()
        for (u in 1:nrow(temp_scen)) {
            if (u == 1) {
                one_scen <- temp_scen[u, ]
                one_scen <- one_scen[, !(names(temp_scen %in% 'time__'))]
                one_scen_sims <- qi_builder(b_sims = param_sims,
                                            newdata = one_scen)
                temp_sims <- rbind(temp_sims, one_scen_sims)
            }
            temp_sims
        }
    })

    return(sims)
}