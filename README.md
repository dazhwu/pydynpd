# pydynpd: Dynamic panel estimation for Difference and System GMM (generalized method-of-moments)
[![DOI](https://zenodo.org/badge/466146436.svg)](https://zenodo.org/badge/latestdoi/466146436)
[![pypi package](https://img.shields.io/pypi/v/pydynpd?style=plastic)](https://pypi.org/project/pydynpd/)

pydynpd is the first python package to implement Difference and System GMM [1][2][3] to estimate dynamic panel data models.

Below is a typical dynamic panel data model:

![y_{it}=\sum^p_{l=1} \alpha_l y_{i,t-l}+\beta x_{i,t}+\gamma s_{i,t}+u_i+\epsilon_{it}](https://latex.codecogs.com/svg.image?y_{it}=\sum^p_{l=1}&space;\alpha_l&space;y_{i,t-l}&plus;\beta&space;x_{i,t}&plus;\gamma&space;s_{i,t}&plus;u_i&plus;\epsilon_{it})
 
In the equation above, x is a predetermined variable that is potentially correlated with past errors, s is a strictly exogenous variable, and u is fixed effect.

## Features supported:
* Differene and System GMM
* One-step and Two-step estimates
* Robust standard errors. For two-step GMM, the calculation suggested by Windmeijer (2005) is used.
* Hansen over-identification test
* Arellano-Bond test for autocorrelation
* Time dummies
* Collapse GMM instruments to limit instrument proliferation


## Installlation:
``` 
pip install pydynpd
``` 

## usage:
``` 
import pandas as pd
from  pydynpd import regression

df = pd.read_csv("data.csv")
command_str='n L(1/2).n w k  | gmm(n, 2 4) gmm(w, 1 3)  iv(k) | timedumm  nolevel'
mydpd = regression.abond(command_str, df, ['id', 'year'])
``` 
### Function abond(command string, data frame, identifier)
where command string consists of two or three parts. The first two parts are required. <br>
* Part one is a list that starts with dependent variable, followed by independent variables. Lag operators can be used in command string. For example, L2.n means to lag variable n two periods. Shortcut L(1/2).n means lags 1 through 2 of variable n, and is equivalent to  L1.n L2.n. 
* Part two desribes how instruments are created. The grammar in this part is similar to that in Stata package XTabond2. More specifically, GMM(varaible list, min_lag max_lag) indicates that lags min_lag through max_lag of each variable included in list of variables are used to generate instruments. For example, GMM(w k, 1 3) means lags 1 through 3 of variables w and k are treated as instruments. On the other hand, IV(variable list) means each variable on variable list is treated as instruments.
* Part three is optional. It includes the following possible options: 
  * onestep: perform one-step estimation rather than the default two-step estimation.
  * nolevel: only perform difference GMM
  * timedumm: include time dummies
  * collapse: collapse instruments



## result:
``` 
Dynamic panel-data estimation, two-step difference GMM
 Group variable: id             Number of obs = 611    
 Time variable: year            Number of groups = 140 
 Number of instruments = 42                            
+-----------+------------+---------------------+------------+-----------+
|     n     |   coef.    | Corrected Std. Err. |     z      |   P>|z|   |
+-----------+------------+---------------------+------------+-----------+
|    L1.n   | 0.2710675  |      0.1382542      | 1.9606462  | 0.0499203 |
|    L2.n   | -0.0233928 |      0.0419665      | -0.5574151 | 0.5772439 |
|     w     | -0.5668527 |      0.2092231      | -2.7093219 | 0.0067421 |
|     k     | 0.3613939  |      0.0662624      | 5.4539824  | 0.0000000 |
| year_1979 | 0.0011898  |      0.0092322      | 0.1288765  | 0.8974554 |
| year_1980 | -0.0316432 |      0.0116155      | -2.7242254 | 0.0064453 |
| year_1981 | -0.0900163 |      0.0206593      | -4.3571693 | 0.0000132 |
| year_1982 | -0.0996210 |      0.0296036      | -3.3651654 | 0.0007650 |
| year_1983 | -0.0693308 |      0.0404276      | -1.7149347 | 0.0863572 |
| year_1984 | -0.0614505 |      0.0475525      | -1.2922666 | 0.1962648 |
+-----------+------------+---------------------+------------+-----------+
Hansen test of overid. restrictions: chi(32) = 32.666 Prob > Chi2 = 0.434
Arellano-Bond test for AR(1) in first differences: z = -1.29 Pr > z =0.198
Arellano-Bond test for AR(2) in first differences: z = -0.31 Pr > z =0.760
``` 
## Similar packages
The objective of the package is similar to the following open-source packages (note: though Stata is a comercial software, xtabond2 is open sourced): <br>
Package | Language
--- | --- 
xtabond2 | Mata （Stata）
plm | R
panelvar | R
pdynmc | R

To compare pydynpd with similar packages, we performed a benchmark test. More specifically, for each package we run 100 times to estimate the same model with the same data. Their estimates and running times (i.e., total running time of 100 tests) are shown in table below. Scripts of this test are included in the "Benchmark" folder. 
### Benchmarks

estimates   | pydynpd | xtabond2 | plm | panelvar 
--- | --- | --- | --- | --- 
L1.n | 0.9453（0.1430）|0.9453（0.1430）| 0.9933 (0.1465)|0.9453（0.1430）
L2.n | -0.0860 (0.1082) |-0.0860 (0.1082)| -0.1640 (0.1071) |-0.0860 (0.1082)
w | -0.4478 (0.1522) | -0.4478 (0.1522)| 0.0594 (0.0284)|-0.4478 (0.1522)
k | 0.1236 (0.0509) |0.1236 (0.0509)|0.1403 (0.0500)|0.1236 (0.0509)|
const | 1.5631 (0.4993) |1.5631 (0.4993)| ... |1.5631 (0.4993)|
--- | --- | --- | --- | --- | --- 
number of instruments derived|51|51|51|51
Hensen Test|96.44|96.44|105.7|96.44
AR(1) Test|-2.35|-2.35|-1.92|...
AR(2) Test|-1.15|-1.15|-0.12|...
running time (secs) | 8.38 | 6.10  | 14.92 |784.9 

As shown in table above, pydynpd produces consistent results compared with xtabond2 and panelvar. plm has different results because it doesn't include intercept in system GMM in its model. As for pdynmc, we first set "include.x" to be TRUE, but it crashed. Then we set "include.x" to be FALSE, it finished calculation. However, it only included lagged dependent variables. Maybe we didn't use pdynmc properly.

Regarding runnint time, in thoery xtabond2 should have a clear advantage because its calculation part was compiled. The result confirms this. However, developed in pure python, pydynpd is not far behind of xtabond2. Moreover, it is significanly faster than the three R packages which are interpreted scripts just like pydynpd.

## References
<a id="1">[1]</a> 
Arellano, M., & Bond, S. (1991). Some tests of specification for panel data: Monte Carlo evidence and an application to employment equations. The review of economic studies, 58(2), 277-297.

<a id="2">[2]</a> 
Arellano, M., & Bover, O. (1995). Another look at the instrumental variable estimation of error-components models. Journal of econometrics, 68(1), 29-51.

<a id="3">[3]</a> 
Blundell, R., & Bond, S. (1998). Initial conditions and moment restrictions in dynamic panel data models. Journal of econometrics, 87(1), 115-143.

<a id="4">[4]</a>
Roodman, D. (2009). How to do xtabond2: An introduction to difference and system GMM in Stata. The stata journal, 9(1), 86-136.

<a id="5">[5]</a> 
Windmeijer, F. (2005). A finite sample correction for the variance of linear efficient two-step GMM estimators. Journal of econometrics, 126(1), 25-51.
