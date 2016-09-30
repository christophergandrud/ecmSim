ecmSim
==============

[![Build Status](https://travis-ci.org/christophergandrud/ecmSim.svg?branch=master)](https://travis-ci.org/christophergandrud/ecmSim)

Early stage development version of a *possible* R package for simulating 
quantities of interest from Error Correction Models (ECM).

Use the `ecm_builder` function to simulated the quantities of interest and the
`ecm_plot` function (NOT COMPLETED to plot the results.

# Example

Imagine we have a two time series `dv` and `iv`. We want to estimate the 
relationship between these two series using an error correction model. 

Our estimation model could look like this (assuming we have already created
the change and lag variables and they are each in their own vectors or equal
length): 




```r
# Estimate error correction model
m1 <- lm(d_dv ~ lag_dv + lag_iv + d_iv)
```

We then create a data frame of fitted values for the baseline scenario to 
simulate.


```r
baseline_scen <- data.frame(lag_dv = mean(lag_dv, na.rm = TRUE),
                            lag_iv = mean(lag_iv, na.rm = TRUE))
```

We also specify the "shock" to `iv`, the effects of which we want to compare to
the baseline. 


```r
iv_shock <- sd(d_iv, na.rm = TRUE)
```

We now have all of the information we need to simulate the effects estimated in 
the ECM over 20 periods:


```r
m1_sims <- ecm_builder(obj = m1, lag_iv = 'lag_iv', d_iv = 'd_iv',
                       iv_shock = sd(d_iv, na.rm = TRUE),
                       baseline_df = baseline_scen, t_extent = 20)
```

```
## lag_dv not supplied. Assuming first column of baseline_df is the lagged dependent variable:
## 
##       lag_dv
```

```r
# Show a sample of the simmulation output
head(m1_sims)
```

```
##   time__   qi_min qi_median   qi_max is_shocked
## 1      1 6.094381  6.116984 6.219926      FALSE
## 2      2 6.780677  6.990080 7.212796      FALSE
## 3      3 7.359268  7.673848 7.973042      FALSE
## 4      4 7.782659  8.195578 8.567472      FALSE
## 5      5 8.093274  8.596204 9.073168      FALSE
## 6      6 8.321153  8.899521 9.477739      FALSE
```

The simulated quantity of interest is the value of `dv` at each time point 
(i.e. lagged `dv` + the change in `dv` from the previous period).

We can plot the results (note, in the future there will be a `ecm_plot` function
to simplify this process):


```r
ggplot(m1_sims, aes(time__, qi_median, group = is_shocked, 
                    colour == is_shocked, fill = is_shocked)) +
    geom_line(aes(color = is_shocked)) +
    geom_ribbon(aes(ymin = qi_min, ymax = qi_max), alpha = 0.2) +
    scale_y_continuous(limits = c(0, 25)) + 
    xlab('\nSimulation Time') + ylab('Predicted dv\n') +
    theme_bw()
```

![plot of chunk non-interactive-plot](figure/non-interactive-plot-1.png)




