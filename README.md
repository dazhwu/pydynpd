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
mydpd = regression.pydynpd('n L(1/2).n w k  | gmm(n, 2 4) gmm(w, 1 3)  iv(k) | timedumm  nolevel', df, ['id', 'year'])
``` 


