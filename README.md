# Competing Risks Metrics

In this repository we have metrics for survival analysis with competing risks; calibration (calibration plots and numerical measures), discrimination (concordance index), prediction error (brier score, weighted brier score, scaled brier score and integrated brier score) and clinical utiliy (decision curves). 

The idea of rewriting these metrics is to make them independent on the prediction model. Packages such riskRegression or pec that compute some of these metrics require models from these packages. 
The consequence is that we cannot evaluate models outside of the packages, limiting evaluation performace of new models not included in these packages. 

In the Rmarkdown files it is demostrated that we get same results than in the packages by using our own custom functions. riskRegression or pec are used for comparison with the custom functions. 

We will follow the paper and its dataset: https://www.bmj.com/content/377/bmj-2021-069249

https://github.com/survival-lumc/ValidationCompRisks/blob/main/Prediction_CSC.md

Needed libraries to run the markdowns

```r
library(survival)
library(riskRegression)
library(prodlim)
library(ggplot2)
library(geepack)
library(dplyr)
library(tidyr)
library(pec)
library(rsample)
```

Needed libraries for custom functions:

```r
library(prodlim)
library(ggplot2)
```

Not tested yet for single risk cases. 

Next steps will be confidence intervals computation (bootstrap and crossvalidation), options for ties in concordance index, maybe allow for informative censoring (with censoring probabilities being computed with covariates and not just Aalen Johansen)
