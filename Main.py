# import vaex as va
import pandas as pd
from  pydynpd import regression



import time

#start = time.time()
df = pd.read_csv("data.csv")

mydpd = regression.abond('n L(1/2).n w k  | gmm(n, 2 4) gmm(w, 1 3)  iv(k) | timedumm  nolevel', df, ['id', 'year'])



#print(time.time()-start)
