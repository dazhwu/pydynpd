# import vaex as va
import pandas as pd
from  pydynpd import regression



#import time

#start = time.time()
#for i in range(100):
df = pd.read_csv("data.csv")

command_str='n L(1/2).n w k  | gmm(n, 2 4) gmm(w, 1 3)  iv(k) '
mydpd = regression.abond(command_str, df, ['id', 'year'])



#print(time.time()-start)
