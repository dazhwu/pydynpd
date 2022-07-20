# import vaex as va
import pandas as pd
from  pydynpd import regression



import time


df = pd.read_csv("test_data.csv")


#command_str='y L(1:?).y L(1:?).x  | gmm(y, 2:3) iv(L(1:1).x)| timedumm'
#mydpd = regression.abond(command_str, df, ['id', 'year'])
df = pd.read_csv("data.csv")
#mydpd = regression.abond('n L(1:2).n w k  | gmm(n, 2:4) gmm(w, 1:3)  iv(k) |nolevel fod ', df, ['id', 'year'])
mydpd = regression.abond('n L(1:2).n w k | gmm(n, 2:4) gmm(w, 1:3) iv(k)', df, ['id', 'year'])





