ecmSim
==============

[![Build Status](https://travis-ci.org/christophergandrud/ecmSim.svg?branch=master)](https://travis-ci.org/christophergandrud/ecmSim)

Early stage development version of a *possible* R package for simulating
quantities of interest from Error Correction Models (ECM). This includes
interactions [!INCOMPLETE!].

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
library(ecmSim)

m1_sims <- ecm_builder(obj = m1, lag_iv = 'lag_iv', d_iv = 'd_iv',
                       iv_shock = sd(d_iv, na.rm = TRUE),
                       baseline_df = baseline_scen, t_extent = 20)

# Show a sample of the simmulation output
head(m1_sims)
```

```
##        qi_ time__    qi_min qi_median   qi_max is_shocked
## 1 1.714924      1 1.4201657 1.6273636 1.825654      FALSE
## 2 1.714924      2 1.4201657 1.6273636 1.825654      FALSE
## 3 1.466063      3 1.1353166 1.3414862 1.550681      FALSE
## 4 1.253315      4 0.8708672 1.1078537 1.380757      FALSE
## 5 1.071441      5 0.6591518 0.9163724 1.233472      FALSE
## 6 0.915959      6 0.4974704 0.7549775 1.116688      FALSE
```

The simulated quantity of interest is the change in the value of `dv`. We could
alternatively supply `qi_d_dv = FALSE` to `ecm_builder` to return the 
dependent variable at each time point (i.e. lagged `dv` + the change in `dv` 
from the previous period).

We can plot the results (note, in the future there will be a `ecm_plot` function
to simplify this process):


```r
library(ggplot2)

ggplot(m1_sims, aes(time__, qi_median, group = is_shocked,
                    colour == is_shocked, fill = is_shocked)) +
    geom_line(aes(color = is_shocked)) +
    geom_ribbon(aes(ymin = qi_min, ymax = qi_max), alpha = 0.2) +
    xlab('\nSimulation Time') + ylab('Predicted dv Change\n') +
    theme_bw()
```

![plot of chunk non-interactive-plot](figure/non-interactive-plot-1.png)


# To-do

- [ ] Create `ecm_plot` function.

- [ ] Interaction example.

# See Also

- See [Warner (2016)](http://static1.squarespace.com/static/5555d102e4b01c8e639df2ca/t/57dcad86d482e9d2d5628489/1474080150275/Warner-Conditional-Relationships.pdf)
for details.

- King, Gary, Michael Tomz, and Jason Wittenberg. 2000. "Making the Most of
Statistical Analyses: Improving Interpretation and Presentation." American
Journal of Political Science 44(2): 341-55.
