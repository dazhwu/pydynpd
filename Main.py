# import vaex as va
import pandas as pd
from regression import pydynpd

import time

start = time.time()
df = pd.read_csv("data.csv")

#mydpd = pydpd('n L(1/2).n L(0/1).w k L(0/1).ys | gmm(n, 2 99) iv(w L1.w k ys L1.ys)| twostep robust', df, ['id', 'year'])
mydpd = pydynpd('n L(1/2).n w k  | gmm(n, 2 4) gmm(w k, 1 3)  | timedumm ', df, ['id', 'year'])

print(time.time() - start)


