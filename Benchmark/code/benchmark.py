import time

import pandas as pd

from pydynpd import regression

start = time.time()
for i in range(100):
    df = pd.read_csv("data.csv")

    mydpd = regression.abond('n L(1:2).n w k  | gmm(n, 2:4) gmm(w, 1:3)  iv(k)', df, ['id', 'year'])

print(time.time() - start)
