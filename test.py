import pandas as pd
from  pydynpd import regression

df = pd.read_csv("tourism_covid_data-total.csv")  #, index_col=False)
df['monthly_cases']=df['monthly cases']
#command_str='tourism_demand L1.tourism_demand L1.monthly_cases  | gmm(tourism_demand, 2 6) iv(L1.monthly_cases)| collapse '
command_str='tourism_demand L1.tourism_demand monthly_cases  | gmm(tourism_demand, 2 6) iv(monthly_cases)| nolevel collapse '
mydpd = regression.abond(command_str, df, ['Country', 'month_year'])