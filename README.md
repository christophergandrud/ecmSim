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
##   time__    qi_min qi_median    qi_max is_shocked
## 1      1  8.135323  8.177056  8.335014      FALSE
## 2      2  9.258256  9.652599  9.995788      FALSE
## 3      3 10.258438 10.857229 11.397112      FALSE
## 4      4 11.020414 11.818322 12.556660      FALSE
## 5      5 11.600914 12.581432 13.513218      FALSE
## 6      6 12.043161 13.181285 14.302321      FALSE
```
