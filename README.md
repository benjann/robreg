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

    23jun2025 (version 2.0.9):
    - in Stata 17 or newer, the -svy- prefix is now also allowed with -robreg s- and
      -robreg mm-
    - -robreg lts- now supports variance estimation/SEs and prediction of influence
      functions

    18sep2021 (version 2.0.8):
    - various changes in code related to absorb()/ivar(); some of them due to
      changes in moremata's mm_areg()
    - -robreg ls, absorb()/ivar()- no longer saves fixed effects in e() by default
      to save computer time and memory; type -usave- to save the fixed effects;
      option -nou- has been renamed to -nousave-

    06sep2021 (version 2.0.7):
    - ivar()/absorb() option (fixed effects) now also supported by -robreg s-
    - ivar() now automatically clusters SEs on the group id
    - some minor buf fixes

    02sep2021 (version 2.0.6):
    - -robreg m- with ivar() or absorb() now uses mm_aqreg() to obtain LAD starting
      values; this ensures regression equivariance (unlike the previously employed
      approach based on group de-medianing)

    01sep2021 (version 2.0.5):
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
