# pydynpd: Dynamic panel estimation for Difference and System GMM (generalized method-of-moments)
[![DOI](https://zenodo.org/badge/466146436.svg)](https://zenodo.org/badge/latestdoi/466146436)
[![pypi package](https://img.shields.io/pypi/v/pydynpd?style=plastic)](https://pypi.org/project/pydynpd/)

pydynpd is the first python package to implement Difference and System GMM [1][2][3] to estimate dynamic panel data models.

Below is a typical dynamic panel data model:

![y_{it}=\sum^p_{l=1} \alpha_l y_{i,t-l}+\beta x_{i,t}+\gamma s_{i,t}+u_i+\epsilon_{it}](https://latex.codecogs.com/svg.image?y_{it}=\sum^p_{l=1}&space;\alpha_l&space;y_{i,t-l}&plus;\beta&space;x_{i,t}&plus;\gamma&space;s_{i,t}&plus;u_i&plus;\epsilon_{it})
 
In the equation above, x is a predetermined variable that is potentially correlated with past errors, s is a strictly exogenous variable, and u is fixed effect.

## Features supported:
* Differene and System GMM
* One-step, two-step, and iterative estimates
* First-difference and forward orthogonal deviation transformations
* Robust standard errors. For two-step GMM, the calculation suggested by Windmeijer (2005) is used.
* Hansen over-identification test
* Arellano-Bond test for autocorrelation
* Time dummies
* Collapse GMM instruments to limit instrument proliferation
* Search for models based on users' request, rather than just run the model specified by users as other packages do


## Installation:
``` 
pip install pydynpd
``` 
This package requires: numpy, scipy, pandas, and PrettyTable

## Usage:
``` 
import pandas as pd
from  pydynpd import regression

df = pd.read_csv("data.csv")
command_str='n L(1:2).n w k  | gmm(n, 2:4) gmm(w, 1:3)  iv(k) | timedumm  nolevel'
mydpd = regression.abond(command_str, df, ['id', 'year'])
``` 
result:
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
## Tutorial
A detailed tutorial is given in the following two documents:<br>
[input of the abond command](https://github.com/dazhwu/pydynpd/blob/main/vignettes/Tutorial.ipynb).<br>
[output of the abond command](https://github.com/dazhwu/pydynpd/blob/main/vignettes/API.md).

## Similar packages
The objective of the package is similar to the following open-source packages: <br>
Package | Language | version
--- | --- | ---
plm | R | 2.6-1
panelvar | R| 0.5.3
pdynmc | R| 0.9.7

To compare pydynpd with similar packages, we performed performance tests. More specifically, in each test for each package we run 100 times to estimate the same model with the same data. For verification, the tests also include Stata package xtabond2 though Stata is a commercial software. We use xtabond2 for regression result verification because it is the most popular package in estimating dynamic panel models. Figure below is from one of the tests. Note that directly comparing xtabond2's speed with R or Python packages is a little unfair because the calculation part of xtabond2 was compiled while pydynpd and the three R packages are interpreted; xtabond2 should have a clear advantage on speed. 

![Alt text](https://raw.githubusercontent.com/dazhwu/pydynpd/main/Benchmark/images/Test_1.svg)

Though developed in pure python, pydynpd is not far behind of xtabond2. Moreover, it is significanly faster than the three R packages which are interpreted scripts just like pydynpd.

A detailed description of the tests can be found [here](https://github.com/dazhwu/pydynpd/blob/main/Benchmark/performance_comparison.md)

## FAQs
### How to extract coefficients from regression?
For example, if you run:
```
df = pd.read_csv("data.csv")
mydpd = regression.abond('n L(1:2).n w k  | gmm(n, 2:4) gmm(w, 1:3)  iv(k)  ', df, ['id', 'year'])
```

The output regression table will be 
```
+------+------------+---------------------+------------+-----------+-----+
|  n   |   coef.    | Corrected Std. Err. |     z      |   P>|z|   |     |
+------+------------+---------------------+------------+-----------+-----+
| L1.n | 0.9453810  |      0.1429764      | 6.6121470  | 0.0000000 | *** |
| L2.n | -0.0860069 |      0.1082318      | -0.7946553 | 0.4268140 |     |
|  w   | -0.4477795 |      0.1521917      | -2.9422068 | 0.0032588 |  ** |
|  k   | 0.1235808  |      0.0508836      | 2.4286941  | 0.0151533 |  *  |
| _con | 1.5630849  |      0.4993484      | 3.1302492  | 0.0017466 |  ** |
+------+------------+---------------------+------------+-----------+-----+
```
If you want to programably extract a value, for example, the first z value (6.6121470) then you can add the following:
```
>>>mydpd.models[0].regression_table.iloc[0]['z_value']
6.6121469997085915
```
Basically, the object mydpd returned above contains models because pydynpd allows us to run and compare multiple models at the same time. By default, it only contains one model which is models[0]. A model has a regression table which is a pandas dataframe:
```
 >>>mydpd.models[0].regression_table

  variable  coefficient   std_err   z_value       p_value  sig
0     L1.n     0.945381  0.142976  6.612147  3.787856e-11  ***
1     L2.n    -0.086007  0.108232 -0.794655  4.268140e-01     
2        w    -0.447780  0.152192 -2.942207  3.258822e-03   **
3        k     0.123581  0.050884  2.428694  1.515331e-02    *
4     _con     1.563085  0.499348  3.130249  1.746581e-03   **

```
So you can extract any value from this dataframe.

### How to use pydynpd with R?
First, you need to install Python on your computer; then install pydynpd.
```
pip install pydynpd
```
Second, in R environment install package reticulate:
```
install.packages("reticulate")
```
Third, you configure Rstudio so that it can communicate with Python installed in step 1. You can find instruction at
https://www.rstudio.com/blog/rstudio-v1-4-preview-python-support/

Finally, you can use the following template to call pydynpd from R. For comparision, the corresponding Python code is also incuded.
<table>
 <tr>
  <td>   R  </td>
  <td>   <pre lang="R">
library(reticulate) 
dynpd <- import("pydynpd.regression", convert = TRUE)
fd <- import("pandas", convert=TRUE)
df <- fd$read_csv("data.csv")

result <- dynpd$abond('n L(1:2).n w k | gmm(n, 2:4) gmm(w, 1:3) iv(k)', df, c('id', 'year'))
</pre>
</td>  </tr>
 <tr> 
  <td>    Python   </td>
  <td>    <pre lang="Python">
import pandas as pd
from  pydynpd import regression
df = pd.read_csv("data.csv")
mydpd = regression.abond('n L(1:2).n w k | gmm(n, 2:4) gmm(w, 1:3) iv(k)', df, ['id', 'year'])
</pre>
  </td>
</tr>
</table>

Code above generates the following result:
```
 Dynamic panel-data estimation, two-step system GMM
 Group variable: id                               Number of obs = 611     
 Time variable: year                              Min obs per group: 4    
 Number of instruments = 51                       Max obs per group: 6    
 Number of groups = 140                           Avg obs per group: 4.36 
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
As you can see, you don't need to change the command string in R. The only parameter you have to change is the identifiers; ['id', 'year'] in Python is changed to c('id', 'year') in R. Also, from R you can access the properties of the result above the same way you work on Python. For example, after running code above if you run the following R script:
```
reg_table=result$models[[1]]$regression_table
print(reg_table)
```
The output is:
```
  variable coefficient    std_err    z_value      p_value sig
1     L1.n  0.94538100 0.14297640  6.6121470 3.787856e-11 ***
2     L2.n -0.08600694 0.10823176 -0.7946553 4.268140e-01    
3        w -0.44777955 0.15219173 -2.9422068 3.258822e-03  **
4        k  0.12358078 0.05088363  2.4286941 1.515331e-02   *
5     _con  1.56308487 0.49934839  3.1302492 1.746581e-03  **
```
In the example above, reg_table is an R data frame.

## Contributing
There are several ways to contribute to pydynpd:

Submit issue/bug reports [here](https://github.com/dazhwu/pydynpd/issues/), or try to fix the problem yourself and then submit a [pull request](https://github.com/dazhwu/pydynpd/pulls).

Browse the source code and see if anything looks out of place - let us know!

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
