# plm
```
library(plm)
abdata=read.csv("data.csv")
pd <- pdata.frame(abdata, index = c("id", "year"), drop.index = TRUE)
z1<-pgmm(n ~ 1+ lag(n, 1:2) + w + k |lag(n, 2:4) + lag(w, 1:3), data=pd, effect='individual',
         model="twosteps" ,transformation='ld', robust = TRUE)
summary(z1)

```

```
pgmm(formula = n ~ lag(n, 1:2) + w + k | lag(n, 2:4) + lag(w, 
    1:3), data = pd, effect = "individual", model = "twosteps", 
    transformation = "ld", robust = TRUE)

Unbalanced Panel: n = 140, T = 7-9, N = 1031

Number of Observations Used: 1362
Residuals:
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
-1.716999 -0.039468  0.000000  0.001151  0.049452  1.057841 

Coefficients:
              Estimate Std. Error z-value  Pr(>|z|)    
lag(n, 1:2)1  0.993296   0.146555  6.7776 1.222e-11 ***
lag(n, 1:2)2 -0.164000   0.107125 -1.5309  0.125791    
w             0.059379   0.028402  2.0906  0.036560 *  
k             0.140340   0.050027  2.8053  0.005027 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Sargan test: chisq(47) = 105.7369 (p-value = 2.0581e-06)
Autocorrelation test (1): normal = -1.926883 (p-value = 0.053994)
Autocorrelation test (2): normal = -0.1281159 (p-value = 0.89806)
Wald test for coefficients: chisq(4) = 8031.159 (p-value = < 2.22e-16)


```

# panelvar
```
library(panelvar)
abdata=read.csv("data.csv")

p1 <-pvargmm(
  dependent_vars = c("n"),
  lags = 2,
  predet_vars = c("w"),
  exog_vars=c("k"),
  transformation = "fd",
  data = abdata,
  panel_identifier = c("id", "year"),
  steps = c("twostep"),
  system_instruments = TRUE,
  max_instr_dependent_vars = 3,
  max_instr_predet_vars = 3,
  min_instr_dependent_vars = 1L,
  min_instr_predet_vars = 1L,
  collapse = FALSE,
  progressbar = FALSE
)
summary(p1)

```

```
-------------------------------------------------
Dynamic Panel VAR estimation, two-step GMM 
---------------------------------------------------
Transformation: First-differences 
Group variable: id 
Time variable: year 
Number of observations = 611 
Number of groups = 140 
Obs per group: min = 4 
               avg = 4.364286 
               max = 6 
Number of instruments = 51 

===================
        n          
-------------------
lag1_n   0.9454 ***
        (0.1430)   
lag2_n  -0.0860    
        (0.1082)   
w       -0.4478 ** 
        (0.1522)   
k        0.1236 *  
        (0.0509)   
const    1.5631 ** 
        (0.4993)   
===================
*** p < 0.001; ** p < 0.01; * p < 0.05

---------------------------------------------------
Instruments for  equation
 Standard
  FD.(k)
 GMM-type
  Dependent vars: L(2, 4)
  Predet vars: L(1, 3)
  Collapse =  FALSE 
---------------------------------------------------

Hansen test of overid. restrictions: chi2(46) = 96.44 Prob > chi2 = 0
(Robust, but weakened by many instruments.)


```

# pdynmc

```
library(pdynmc)
abdata=read.csv("data.csv")
mc_1 <- pdynmc(dat=abdata,varname.i = "id", varname.t = "year",
               use.mc.diff = TRUE, use.mc.lev = TRUE, use.mc.nonlin = FALSE,
               include.y = TRUE, varname.y = "n", lagTerms.y = 2, maxLags.y=4,
               inst.stata = TRUE, include.x = TRUE,               
               varname.reg.pre = c("w"), lagTerms.reg.pre = c(0), maxLags.reg.pre = c(3),
               fur.con = TRUE, fur.con.diff = TRUE, fur.con.lev = TRUE,
               varname.reg.fur = c("k"),lagTerms.reg.fur = c(0),
               w.mat = "iid.err", std.err = "corrected", estimation = "twostep",
               opt.meth = "none")
summary(mc_1)
mtest.fct(mc_1, order = 2)
```
```
Error in mapply(ti = ti.temp, t.end = tend.temp, lagTerms = lagTerms, : non-numeric argument to binary operator
Traceback:

1. pdynmc(dat = abdata, varname.i = "id", varname.t = "year", use.mc.diff = TRUE, 
 .     use.mc.lev = TRUE, use.mc.nonlin = FALSE, include.y = TRUE, 
 .     varname.y = "n", lagTerms.y = 2, maxLags.y = 4, inst.stata = TRUE, 
 .     include.x = TRUE, varname.reg.pre = c("w"), lagTerms.reg.pre = c(0), 
 .     maxLags.reg.pre = c(3), fur.con = TRUE, fur.con.diff = TRUE, 
 .     fur.con.lev = TRUE, varname.reg.fur = c("k"), lagTerms.reg.fur = c(0), 
 .     w.mat = "iid.err", std.err = "corrected", estimation = "twostep", 
 .     opt.meth = "none")
2. lapply(X = i_cases, FUN = Z_i.fct, Time = Time, varname.i = varname.i, 
 .     use.mc.diff = use.mc.diff, use.mc.lev = use.mc.lev, use.mc.nonlin = use.mc.nonlin, 
 .     use.mc.nonlinAS = use.mc.nonlinAS, include.y = include.y, 
 .     varname.y = varname.y, inst.stata = inst.stata, include.dum = include.dum, 
 .     dum.diff = dum.diff, dum.lev = dum.lev, colnames.dum = colnames.dum, 
 .     fur.con = fur.con, fur.con.diff = fur.con.diff, fur.con.lev = fur.con.lev, 
 .     varname.reg.estParam.fur = varname.reg.estParam.fur, include.x = include.x, 
 .     end.reg = end.reg, varname.reg.end = varname.reg.end, pre.reg = pre.reg, 
 .     varname.reg.pre = varname.reg.pre, ex.reg = ex.reg, varname.reg.ex = varname.reg.ex, 
 .     maxLags.y = maxLags.y, lagTerms.y = lagTerms.y, max.lagTerms = max.lagTerms, 
 .     maxLags.reg.end = maxLags.reg.end, maxLags.reg.pre = maxLags.reg.pre, 
 .     maxLags.reg.ex = maxLags.reg.ex, inst.reg.ex.expand = inst.reg.ex.expand, 
 .     dat = dat, dat.na = dat.na)
3. FUN(X[[i]], ...)
4. do.call(what = "cbind", args = sapply(FUN = LEV.pre.fct, i = i, 
 .     varname.ex.pre.temp, T.mcLev = T.mcLev.temp, use.mc.diff = use.mc.diff, 
 .     inst.stata = inst.stata, Time = Time, varname.i = varname.i, 
 .     lagTerms = max.lagTerms, dat = dat, dat.na = dat.na))
5. sapply(FUN = LEV.pre.fct, i = i, varname.ex.pre.temp, T.mcLev = T.mcLev.temp, 
 .     use.mc.diff = use.mc.diff, inst.stata = inst.stata, Time = Time, 
 .     varname.i = varname.i, lagTerms = max.lagTerms, dat = dat, 
 .     dat.na = dat.na)
6. lapply(X = X, FUN = FUN, ...)
7. FUN(X[[i]], ...)
8. Matrix::bdiag(do.call(what = diag, args = list(mapply(ti = ti.temp, 
 .     t.end = tend.temp, lagTerms = lagTerms, FUN = datLEV.pre.fct, 
 .     varname = varname, MoreArgs = list(i = i, use.mc.diff = use.mc.diff, 
 .         inst.stata = inst.stata, dat = dat, dat.na = dat.na, 
 .         varname.i = varname.i, Time = Time)) * as.vector(!is.na(diff(dat.na[dat.na[, 
 .     varname.i] == i, varname][(lagTerms - 1):Time]))))))
9. do.call(what = diag, args = list(mapply(ti = ti.temp, t.end = tend.temp, 
 .     lagTerms = lagTerms, FUN = datLEV.pre.fct, varname = varname, 
 .     MoreArgs = list(i = i, use.mc.diff = use.mc.diff, inst.stata = inst.stata, 
 .         dat = dat, dat.na = dat.na, varname.i = varname.i, Time = Time)) * 
 .     as.vector(!is.na(diff(dat.na[dat.na[, varname.i] == i, varname][(lagTerms - 
 .         1):Time])))))

```
# pydynpd
```
import pandas as pd
from  pydynpd import regression
df=pd.read_csv("data.csv")

mydpd = regression.abond('n L(1:2).n w k | gmm(n, 2:4) gmm(w, 1:3) iv(k)', df, ['id', 'year'])
```
```
Dynamic panel-data estimation, two-step system GMM
 Group variable: id             Number of obs = 751     
 Time variable: year            Min obs per group: 5    
 Number of instruments = 51     Max obs per group: 7    
 Number of groups = 140         Avg obs per group: 5.36 
+------+------------+---------------------+------------+-----------+-----+
|  n   |   coef.    | Corrected Std. Err. |     z      |   P>|z|   |     |
+------+------------+---------------------+------------+-----------+-----+
| L1.n | 0.9453810  |      0.1429764      | 6.6121470  | 0.0000000 | *** |
| L2.n | -0.0860069 |      0.1082318      | -0.7946553 | 0.4268140 |     |
|  w   | -0.4477795 |      0.1521917      | -2.9422068 | 0.0032588 |  ** |
|  k   | 0.1235808  |      0.0508836      | 2.4286941  | 0.0151533 |  *  |
| _con | 1.5630849  |      0.4993484      | 3.1302492  | 0.0017466 |  ** |
+------+------------+---------------------+------------+-----------+-----+
Hansen test of overid. restrictions: chi(46) = 96.442 Prob > Chi2 = 0.000
Arellano-Bond test for AR(1) in first differences: z = -2.35 Pr > z =0.019
Arellano-Bond test for AR(2) in first differences: z = -1.15 Pr > z =0.251
```



# xtabond2 (default)
```
insheet using "data.csv"
xtset(id year)
xtabond2 n L(1/2).n w k , gmm(n, lag(2 4)) gmm(w, lag(1 3)) iv(k ) nolevel twostep robust 
```

```
Favoring space over speed. To switch, type or click on mata: mata set matafavor speed, perm.
Warning: Two-step estimated covariance matrix of moments is singular.
  Using a generalized inverse to calculate optimal weighting matrix for two-step estimation.
  Difference-in-Sargan/Hansen statistics may be negative.

Dynamic panel-data estimation, two-step difference GMM
------------------------------------------------------------------------------
Group variable: id                              Number of obs      =       611
Time variable : year                            Number of groups   =       140
Number of instruments = 36                      Obs per group: min =         4
Wald chi2(0)  =         .                                      avg =      4.36
Prob > chi2   =         .                                      max =         6
------------------------------------------------------------------------------
             |              Corrected
           n |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
           n |
         L1. |   .1700616   .1046652     1.62   0.104    -.0350784    .3752016
         L2. |  -.0113381   .0377205    -0.30   0.764    -.0852688    .0625926
             |
           w |  -.9510582   .1277298    -7.45   0.000    -1.201404   -.7007124
           k |   .4637223   .0718328     6.46   0.000     .3229325    .6045121
------------------------------------------------------------------------------
Instruments for first differences equation
  Standard
    D.k
  GMM-type (missing=0, separate instruments for each period unless collapsed)
    L(1/3).w
    L(2/4).n
------------------------------------------------------------------------------
Arellano-Bond test for AR(1) in first differences: z =  -1.19  Pr > z =  0.235
Arellano-Bond test for AR(2) in first differences: z =  -0.81  Pr > z =  0.417
------------------------------------------------------------------------------
Sargan test of overid. restrictions: chi2(32)   =  91.61  Prob > chi2 =  0.000
  (Not robust, but not weakened by many instruments.)
Hansen test of overid. restrictions: chi2(32)   =  47.86  Prob > chi2 =  0.035
  (Robust, but weakened by many instruments.)

Difference-in-Hansen tests of exogeneity of instrument subsets:
  gmm(n, lag(2 4))
    Hansen test excluding group:     chi2(15)   =  23.75  Prob > chi2 =  0.069
    Difference (null H = exogenous): chi2(17)   =  24.11  Prob > chi2 =  0.117
  gmm(w, lag(1 3))
    Hansen test excluding group:     chi2(14)   =  17.25  Prob > chi2 =  0.243
    Difference (null H = exogenous): chi2(18)   =  30.61  Prob > chi2 =  0.032
  iv(k)
    Hansen test excluding group:     chi2(31)   =  38.33  Prob > chi2 =  0.171
    Difference (null H = exogenous): chi2(1)    =   9.53  Prob > chi2 =  0.002
```
# xtabond2 (speed)

```
mata: mata set matafavor speed, perm
insheet using "data.csv"
xtset(id year)
xtabond2 n L(1/2).n w k , gmm(n, lag(2 4)) gmm(w, lag(1 3)) iv(k ) nolevel twostep robust 

```
```

Favoring speed over space. To switch, type or click on mata: mata set matafavor space, perm.
Warning: Two-step estimated covariance matrix of moments is singular.
  Using a generalized inverse to calculate optimal weighting matrix for two-step estimation.
  Difference-in-Sargan/Hansen statistics may be negative.

Dynamic panel-data estimation, two-step difference GMM
------------------------------------------------------------------------------
Group variable: id                              Number of obs      =       611
Time variable : year                            Number of groups   =       140
Number of instruments = 36                      Obs per group: min =         4
Wald chi2(0)  =         .                                      avg =      4.36
Prob > chi2   =         .                                      max =         6
------------------------------------------------------------------------------
             |              Corrected
           n |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
           n |
         L1. |   .1700616   .1046652     1.62   0.104    -.0350784    .3752016
         L2. |  -.0113381   .0377205    -0.30   0.764    -.0852688    .0625926
             |
           w |  -.9510582   .1277298    -7.45   0.000    -1.201404   -.7007124
           k |   .4637223   .0718328     6.46   0.000     .3229325    .6045121
------------------------------------------------------------------------------
Instruments for first differences equation
  Standard
    D.k
  GMM-type (missing=0, separate instruments for each period unless collapsed)
    L(1/3).w
    L(2/4).n
------------------------------------------------------------------------------
Arellano-Bond test for AR(1) in first differences: z =  -1.19  Pr > z =  0.235
Arellano-Bond test for AR(2) in first differences: z =  -0.81  Pr > z =  0.417
------------------------------------------------------------------------------
Sargan test of overid. restrictions: chi2(32)   =  91.61  Prob > chi2 =  0.000
  (Not robust, but not weakened by many instruments.)
Hansen test of overid. restrictions: chi2(32)   =  47.86  Prob > chi2 =  0.035
  (Robust, but weakened by many instruments.)

Difference-in-Hansen tests of exogeneity of instrument subsets:
  gmm(n, lag(2 4))
    Hansen test excluding group:     chi2(15)   =  23.75  Prob > chi2 =  0.069
    Difference (null H = exogenous): chi2(17)   =  24.11  Prob > chi2 =  0.117
  gmm(w, lag(1 3))
    Hansen test excluding group:     chi2(14)   =  17.25  Prob > chi2 =  0.243
    Difference (null H = exogenous): chi2(18)   =  30.61  Prob > chi2 =  0.032
  iv(k)
    Hansen test excluding group:     chi2(31)   =  38.33  Prob > chi2 =  0.171
    Difference (null H = exogenous): chi2(1)    =   9.53  Prob > chi2 =  0.002

```

