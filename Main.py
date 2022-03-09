# import vaex as va
import pandas as pd
from dynamic_panel_data_model import dpd


import time



start = time.time()
df = pd.read_csv("data.csv")

#command=command('n L(1/2).n w k| gmm(n, 2 4) gmm(w, 1 3) iv(k)| twostep nolevel robust')
#variables=command.variables

mydpd = dpd('n L(1/2).n w k| gmm(n, 2 4) gmm(w, 1 3) iv(k)| twostep nolevel robust', df, ['id', 'year'])

#data=panel_data(df, ['id','year'], variables, balanced=True)

#prepare_reg()

#mydpd = dpd('n L(1/2).n w k| gmm(n, 2 4) gmm(w, 1 3) iv(k)| twostep nolevel robust', df, ['id', 'year'])


print(time.time() - start)


