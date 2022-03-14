# pydynpd: Dynamic panel estimation for Difference and System GMM (generalized method-of-moments)
[![DOI](https://zenodo.org/badge/466146436.svg)](https://zenodo.org/badge/latestdoi/466146436)
[![pypi package](https://img.shields.io/pypi/v/pydynpd?style=plastic)](https://pypi.org/project/pydynpd/)

pydynpd is the first python package to implement Difference and System GMM [1][2][3] to estimate dynamic panel data models.

Below is a typical dynamic panel data model:

![y_{it}=\sum^p_{l=1} \alpha_l y_{i,t-l}+\beta x_{i,t}+\gamma s_{i,t}+u_i+\epsilon_{it}](https://latex.codecogs.com/svg.image?y_{it}=\sum^p_{l=1}&space;\alpha_l&space;y_{i,t-l}&plus;\beta&space;x_{i,t}&plus;\gamma&space;s_{i,t}&plus;u_i&plus;\epsilon_{it})
 
In the equation above, x is a predetermined variable that is potentially correlated with past errors, s is a strictly exogenous variable, and u is fixed effect.

## Features supported:
* Differene and System GMM
* One-step and Two-step
* Robust standard errors. For two-step GMM, the calculation suggested by Windmeijer (2005) is used.
* Hansen over-identification test
* Arellano-Bond test for autocorrelation
* Time dummies


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
## Benchmarks
The objective of the package is similar to the following packages: <br>
Package | Language
--- | --- 
xtabond2 | Stata
plm | R
panelvar | R
pdynmc | R

Stata command reghdfe and the R function felm. The package tends to be much faster than these two options.

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
