# pydynpd: Dynamic panel estimation for Difference and System GMM (generalized method-of-moments)
[![DOI](https://zenodo.org/badge/466146436.svg)](https://zenodo.org/badge/latestdoi/466146436)
[![pypi package](https://img.shields.io/pypi/v/pydynpd?style=plastic)](https://pypi.org/project/pydynpd/)


Installlation: <br>
pip install pydynpd <br>

usage: <br>
``` 
import pandas as pd
from  pydynpd import regression

df = pd.read_csv("data.csv")
mydpd = regression.abond('n L(1/2).n w k  | gmm(n, 2 4) gmm(w, 1 3)  iv(k) | timedumm  nolevel', df, ['id', 'year'])
``` 

result: <br>
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
