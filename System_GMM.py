import numpy as np
from pandas import DataFrame

import panel_data as diff_gmm


class system_gmm():

    def __init__(self, df: DataFrame, identifiers, variables):
        diff_result = diff_gmm.panel_data(df, ['id', 'year'], variables)
        print(len(diff_result.z_list))

