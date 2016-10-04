#' Create simulations for long-term effects from error correction models (ECM)
#'
#' @param obj a fitted model object from an ECM model estimated with
#' \code{\link{lm}}.
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
#' @param lag_iv_2 character string identifying variable name of the second
#' lagged independent variable in an interaction with \code{lag_iv}, if
#' applicable.
#' @param d_iv_2 character string identifying variable name of the second
#' shocked independent variable in an interaction with \code{lag_iv}, if
#' applicable.
#' @param lag_iv_lag_iv2_interaction character string identifying the interaction
#' term for \code{lag_iv} * \code{lag_iv_2}, if applicable.
#' @param d_iv_d_iv2_interaction character string identifying the interaction
#' term for \code{d_iv} * \code{d_iv_2}, if applicable.
#' @param lag_iv_d_iv2_interaction character string identifying the interaction
#' term for \code{lag_iv} * \code{d_iv_2}, if applicable.
#' @param d_iv_lag_iv2_interaction character string identifying the interaction
#' term for \code{d_iv} * \code{lag_iv_2}, if applicable.
#' @param iv_2_shock numeric shock to \code{iv_2}. If not specified, set to
#' 0.
#' @param t_extent numeric specifying the time points from the shock to
#' simulate the long-term effects of the shock for. The default is 5 time
#' points.
#' @param nsim numeric. Number of simulations to draw.
#' @param ci the proportion of the central interval of the simulations to
#' return. Must be in (0, 1] or equivalently (0, 100]. Note: if \code{ci = 1}
#' then the full interval (i.e. 100 percent) is assumed.
#' @param slim logical indicating whether to (if \code{FALSE}) return all
#' simulations in the central interval specified by \code{ci} for each fitted
#' scenario or (if \code{TRUE}) just the minimum, median, and maxium values.
#' See \code{\link{qi_slimmer}} for more details.
#' @param mu an optional vector giving the means of the variables estimated
#' from an ECM model. If \code{obj} is supplied then \code{mu} is ignored.
#' @param Sigma an optional positive-definite symmetric matrix estimated
#' from an ECM model. The matrix specifies the covariance matrix of the
#' variables. If \code{obj}  is supplied then \code{Sigma} is ignored.
#' Note that if the model includes an intercept the corresponding column in
#' \code{Sigma} must be called \code{intercept_}.
#'
#' @details If an interaction with the shocked variable is included, the
#' function assumes a "general model" setup as discussed in Warner (2016).
#'
#'
#' @examples
#' # Simulated data to estimate the model from
#' iv <- rnorm(10)
#' lag_iv <- c(NA, iv[-length(iv)])
#' d_iv <-  iv - lag_iv
#'
#' b0 = 2
#' b1 = -0.2
#' b2 = 1.1
#' b3 = 2.2
#' starting_d_dv <- -4.1
#' eps <- rnorm(iv, sd = 0.2)
#'
#' dv <- vector()
#' d_dv <- vector()
#' lag_dv <- vector()
#' for (i in 2:length(lag_iv)) {
#'     if (i == 2) {
#'         dv_start <- b0 + b2*lag_iv[i] + b3*d_iv[i] + eps[i]
#'                           lag_dv[3] <- dv_start + starting_d_dv
#'     }
#'     else if (i == 3) {
#'         d_dv[i] <- b0 + b1*lag_dv[i] + b2*lag_iv[i] + b3*d_iv[i] + eps[i]
#'     }
#'     else if (i > 3) {
#'         lag_dv[i] <- lag_dv[i-1] + d_dv[i-1]
#'         d_dv[i] <- b0 + b1*lag_dv[i] + b2*lag_iv[i] + b3*d_iv[i] + eps[i]
#'     }
#' }
#'
#' # Estimate error correction model
#' m1 <- lm(d_dv ~ lag_dv + lag_iv + d_iv)
#'
#' # Create baseline scenario
#' baseline_scen <- data.frame(lag_dv = mean(lag_dv, na.rm = TRUE),
#'                             lag_iv = mean(lag_iv, na.rm = TRUE))
#'
#' # Simulate effects of a 1 standard deviation shock to the IV over 20 time
#' ## periods
#' m1_sims <- ecm_builder(obj = m1, lag_iv = 'lag_iv', d_iv = 'd_iv',
#'                        iv_shock = sd(d_iv, na.rm = TRUE),
#'                        baseline_df = baseline_scen, t_extent = 20)
#'
#'
#' @importFrom coreSim b_sim qi_builder qi_slimmer
#' @export

ecm_builder <- function(obj, baseline_df, lag_dv,
                        lag_iv, d_iv, iv_shock,
                        lag_iv_2, d_iv_2,
                        lag_iv_lag_iv2_interaction, d_iv_d_iv2_interaction,
                        lag_iv_d_iv2_interaction, d_iv_lag_iv2_interaction,
                        iv_2_shock,
                        t_extent = 5,
                        nsim = 1000, ci = 0.95, slim = TRUE,
                        mu, Sigma)
{
    time__ <- NULL

    ci <- ecmSim:::ci_check(ci)

    if (t_extent < 3) {
        message('t_extent must be 3 or more time points.\nForcing t_extent = 3.')
        t_extent = 3
    }

    if (missing(lag_dv)) {
        lag_dv <- names(baseline_df)[1]
        message(paste(
            'lag_dv not supplied. Assuming first column of baseline_df is the lagged dependent variable:\n\n     ',
            lag_dv, '\n'))
    }

    baseline_scenario <- ecmSim:::df_repeat(baseline_df, n = t_extent)
    baseline_scenario$time__ <- 1:t_extent
    baseline_scenario[, lag_dv][baseline_scenario$time__ > 1] <- NA

    if (!missing(lag_iv_2) & !missing(lag_iv_lag_iv2_interaction))
        baseline_scenario[,
            lag_iv_lag_iv2_interaction] <- baseline_scenario[, lag_iv] *
                                            baseline_scenario[, lag_iv_2]

    non_lag_dv_names <- names(baseline_scenario)[!(names(baseline_scenario) %in%
                                                       c(lag_dv, 'time__'))]

    baseline_scenario$is_shocked <- FALSE

    # Create shock fitted values
    shocked <- ecmSim:::df_repeat(baseline_df, n = t_extent)
    shocked$time__ <- 1:t_extent
    if (missing(iv_shock)) iv_shock <- 0
    shocked[2, d_iv] <- iv_shock
    shocked[is.na(shocked)] <- 0
    shocked[3:t_extent, lag_iv] <- shocked[2, lag_iv] + iv_shock
    shocked[2:t_extent, lag_dv] <- NA

    if (!missing(lag_iv_lag_iv2_interaction) &
        !missing(d_iv_d_iv2_interaction) &
        !missing(lag_iv_2) & !missing(d_iv_2) &
        !missing(lag_iv_d_iv2_interaction) &
        !missing(d_iv_lag_iv2_interaction))
    {
        message('Creating shock interactions...')
        if (missing(iv_2_shock)) iv_2_shock <- 0
        # First period no shock to second variable
        shocked[, d_iv_2] <- 0

        # Second period+ shock
        shocked[2, d_iv_2] <- iv_2_shock
        shocked[, d_iv_d_iv2_interaction] <- shocked[, d_iv] *
            shocked[, d_iv_2]

        shocked[3:t_extent, lag_iv_2] <- shocked[2, lag_iv_2] + iv_2_shock
        shocked[, lag_iv_lag_iv2_interaction] <- shocked[, lag_iv] *
            shocked[, lag_iv_2]


        shocked[, lag_iv_d_iv2_interaction] <- shocked[, lag_iv] *
                shocked[, d_iv_2]
        shocked[, d_iv_lag_iv2_interaction] <- shocked[, d_iv] *
                shocked[, lag_iv_2]
    }

    non_lag_dv_names_shocked <- names(shocked)[!(names(shocked) %in%
                                                     c(lag_dv, 'time__'))]

    shocked$is_shocked <- TRUE

    scenarios <- list(baseline = baseline_scenario, shocked = shocked)

    # Simulate parameters
    if (!missing(obj))
        param_sims <- b_sim(obj = obj, nsim = nsim)
    else if (!missing(mu) & !missing(Sigma))
        param_sims <- b_sim(mu = mu, Sigma = Sigma, nsim = nsim)

    sims <- list()
    sims <- lapply(seq_along(scenarios), function(x) {
        temp_scen <- scenarios[[x]]
        is_shocked <- any(temp_scen$is_shocked)
        if (is_shocked) col_names = non_lag_dv_names_shocked
        else col_names = non_lag_dv_names

        temp_scen <- ecmSim:::drop_col(temp_scen, 'is_shocked')

        temp_sims <- data.frame()
        for (u in 1:nrow(temp_scen)) {
            one_scen <- temp_scen[u, ]
            one_scen <- one_scen[, !(names(temp_scen) %in% 'time__')]
            if (u == 1) {
                one_scen_sims <- qi_builder(b_sims = param_sims,
                                            newdata = one_scen,
                                            ci = 1, verbose = FALSE)
                one_scen_sims[, lag_dv] <- one_scen_sims[, lag_dv] +
                                            one_scen_sims[, 'qi_']
                one_scen_sims <- ecmSim:::drop_col(one_scen_sims, 'qi_')
                one_scen_sims$time__ <- u
                one_scen_sims$sim_id__ <- 1:nrow(one_scen_sims)
                temp_sims <- rbind(temp_sims, one_scen_sims)
            }
            else if (u > 1) {
                lag_dv_sim_values <- subset(temp_sims, time__ == u - 1)
                lag_dv_sim_values <- data.frame(lag_dv_sim_values[, lag_dv])
                names(lag_dv_sim_values) <- lag_dv
                one_scen_x <- data.frame(one_scen[, col_names])
                colnames(one_scen_x) <- col_names
                temp_scen_updated <- ecmSim:::expand_dfs(lag_dv_sim_values,
                                                         one_scen_x)
                one_scen_sims <- param_sims[, names(temp_scen_updated)] *
                                    temp_scen_updated
                temp_scen_updated[, lag_dv] <- temp_scen_updated[, lag_dv] +
                                    (rowSums(one_scen_sims) +
                                    param_sims[, 'intercept_'])
                temp_scen_updated$time__ <- u
                temp_scen_updated$sim_id__ <- 1:nrow(temp_scen_updated)
                temp_sims <- rbind(temp_sims, temp_scen_updated)
            }
        }
        temp_sims <- temp_sims[, c(lag_dv, 'time__', 'sim_id__')]
        temp_sims$is_shocked <- is_shocked
        temp_sims
    })

    sims <- data.frame(bind_rows(sims))
    sims$scenario_ <- interaction(sims[, c('time__', 'is_shocked')])
    sims <- qi_central_interval(sims, scenario_var = 'scenario_',
                                qi_var = lag_dv, ci = ci)
    if (slim) {
        sims <- qi_slimmer(sims, qi_var = lag_dv)
        sims$is_shocked <- c(baseline_scenario$is_shocked, shocked$is_shocked)
        sims <- drop_col(sims, c('scenario_', lag_dv, 'sim_id__'))
    }
    else sims <- drop_col(sims, 'scenario_')

    return(sims)
}
