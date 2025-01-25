
import pandas as pd
from  pydynpd import regression

import time

a=time.time()


df = pd.read_csv("data.csv")
for i in range(0,101):
#mydpd = regression.abond('n L(1:2).n w k | gmm(n, 2:4) gmm(w, 1:3) iv(k)', df, ['id', 'year'])
    mydpd = regression.abond('n L(1:2).n w k | gmm(n, 2:.) pred(w k)', df, ['id', 'year'])

print(time.time()-a)