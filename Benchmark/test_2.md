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
               use.mc.diff = TRUE, use.mc.lev = FALSE, use.mc.nonlin = FALSE,
               include.y = TRUE, varname.y = "n", lagTerms.y = 2, maxLags.y=4,
               inst.stata = TRUE, include.x = TRUE,               
               varname.reg.pre = c("w"), lagTerms.reg.pre = c(0), maxLags.reg.pre = c(3),
               fur.con = TRUE, fur.con.diff = TRUE, fur.con.lev = FALSE,
               varname.reg.fur = c("k"),lagTerms.reg.fur = c(0),
               w.mat = "iid.err", std.err = "corrected", estimation = "twostep",
               opt.meth = "none")
summary(mc_1)
mtest.fct(mc_1, order = 2)
```
```
Dynamic linear panel estimation (twostep)
Estimation steps: 2

Coefficients:
     Estimate Std.Err.rob z-value.rob Pr(>|z.rob|)    
L1.n  0.17078     0.10597       1.611        0.107    
L2.n -0.01186     0.03862      -0.307        0.759    
L0.w -0.96426     0.12689      -7.599       <2e-16 ***
L0.k  0.46357     0.07237       6.406       <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

 36 total instruments are employed to estimate 4 parameters
 35 linear (DIF) 
 1 further controls (DIF) 
 no time dummies 
 
J-Test (overid restrictions):  47.49 with 32 DF, pvalue: 0.0383
F-Statistic (slope coeff):  408.98 with 4 DF, pvalue: <0.001
F-Statistic (time dummies):  no time dummies included in estimation

	Arellano and Bond (1991) serial correlation test of degree 2

data:  2step GMM Estimation
normal = -0.9218, p-value = 0.3566
alternative hypothesis: serial correlation of order 2 in the error terms



	Arellano and Bond (1991) serial correlation test of degree 2

data:  2step GMM Estimation
normal = -0.9218, p-value = 0.3566
alternative hypothesis: serial correlation of order 2 in the error terms

```
# pydynpd
```
import pandas as pd
from  pydynpd import regression
df = pd.read_csv("data.csv")
mydpd = regression.abond('n L(1:2).n w k  | gmm(n, 2:4) gmm(w, 1:3)  iv(k) | nolevel', df, ['id', 'year'])
```
```
Dynamic panel-data estimation, two-step difference GMM
 Group variable: id             Number of obs = 611     
 Time variable: year            Min obs per group: 5    
 Number of instruments = 36     Max obs per group: 7    
 Number of groups = 140         Avg obs per group: 5.36 
+------+------------+---------------------+------------+-----------+-----+
|  n   |   coef.    | Corrected Std. Err. |     z      |   P>|z|   |     |
+------+------------+---------------------+------------+-----------+-----+
| L1.n | 0.1700616  |      0.1046652      | 1.6248154  | 0.1042019 |     |
| L2.n | -0.0113381 |      0.0377205      | -0.3005824 | 0.7637329 |     |
|  w   | -0.9510582 |      0.1277298      | -7.4458585 | 0.0000000 | *** |
|  k   | 0.4637223  |      0.0718328      | 6.4555747  | 0.0000000 | *** |
+------+------------+---------------------+------------+-----------+-----+
Hansen test of overid. restrictions: chi(32) = 47.860 Prob > Chi2 = 0.035
Arellano-Bond test for AR(1) in first differences: z = -1.19 Pr > z =0.235
Arellano-Bond test for AR(2) in first differences: z = -0.81 Pr > z =0.417
```


    command_str='y L1.y L1.x  | gmm(y, 2:4) iv(L1.x)| timedumm '
    mydpd = regression.abond(command_str, df, ['id', 'year'])
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

