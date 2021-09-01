# robreg
Stata module providing robust regression estimators

`robreg` provides a number of robust estimators for linear regression models.
Among them are the high breakdown-point and high efficiency MM estimator, the
Huber and bisquare M estimator, the S estimator, as well as quantile
regression, each supporting robust standard errors based on influence
functions. Furthermore, basic implementations of LMS/LQS (least median/quantile
of squares) and LTS (least trimmed squares) are provided.

Requires: Stata 11 or newer, package `moremata`

To install `robreg` from the SSC Archive, type

    . ssc install robreg, replace
    . ssc install moremata, replace

in Stata.

---

Installation from GitHub:

    . net install robreg, replace from(https://raw.githubusercontent.com/benjann/robreg/main/)
    . net install moremata, replace from(https://raw.githubusercontent.com/benjann/moremata/master/)

---

Main changes:

    25aug2021 (version 2.0.5):
    - new ivar()/absorb() option (fixed effects) in -robreg ls- and -robreg m-
    - -predict ... if 0, ifs- failed; this is fixed
    - fixed alignment of table header in Stata 17

    21apr2021 (version 2.0.4):
    - -robreg hausman- failed if applied to models that were estimated with
      weights; this is fixed
    
    19apr2021 (version 2.0.3):
    - time-series operators are now allowed
    
    16apr2021 (version 2.0.2):
    - options -noheader- and -notable- added
    - option -nodetail- in -robreg hausman- no longer documented; use -notable-
    - robreg can now display results from -xtrobreg-
    
    08apr2021 (version 2.0.1):
    - robreg lts/lqs/lms now store the h-quantile of squared residuals in e(q_h)
      (for use by -predict, subset-)
    
    07apr2021 (version 2.0.0):
    - robreg released on GitHub
