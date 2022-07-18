import math

import numpy as np
from pandas import DataFrame

from pydynpd.common_functions import get_first_diff_table, get_fod_table
from pydynpd.info import options_info


class panel_data():
    def __init__(self, df: DataFrame, identifiers, variables, options: options_info):

        self._individual = identifiers[0]
        self._time = identifiers[1]

        cols = []

        temp_list = [var.name for var in variables['dep_indep'] + variables['iv'] + variables['Dgmm']]
        for var_name in temp_list:
            if var_name not in cols:
                cols.append(var_name)

        # temp_df = df[ [self._individual, self._time] + cols].copy()
        self.N, self.T, self.ids = self.xtset(df, self._individual, self._time)

        if options.timedumm:
            self.col_timedumm = self.add_time_dummy(df, variables, self._time)
        else:
            self.col_timedumm = []

        self.cols = self.ids + cols  # make sure ids is the first column

        self.data = self.make_balanced(df[self.cols + self.col_timedumm].to_numpy(), self.N, self.T)
        num_cols = self.data.shape[1]

        self.fd_data = get_first_diff_table(self.data[:, range(0, num_cols)], self.N)
        if (options.transformation=='fod') & (options.level==False):
            self.fod_data=get_fod_table(self.data, self.N)

    def xtset(self, df: DataFrame, _individual, _time):
        df.sort_values(by=[_individual, _time])
        df['_individual'] = df[_individual].astype('category').cat.codes
        df['_individual'] = df['_individual'].astype('int64')
        N = df['_individual'].unique().size

        df['_time'] = df[_time].astype('category').cat.codes
        df['_time'] = df['_time'].astype('int64')
        T = df['_time'].unique().size

        df['_NT'] = df['_individual'] * T + df['_time']

        if N <= T:
            print(
                'Warning: system and difference GMMs do not work well on long (T>=N) panel data')
        return (N, T, ['_NT'])

    def make_balanced(self, ori, n_individual, n_time):
        arr_full = np.empty((n_individual * n_time, ori.shape[1]), dtype='float64')

        arr_full[:] = np.NaN
        # arr_full[:, 0] = np.repeat(range(0, N), T)
        arr_full[:, 0] = range(0, n_individual * n_time)
        # arr_full[:, 1] = arr_full[:, 2] % T

        mask = np.in1d(arr_full[:, 0], ori[:, 0])

        arr_full[mask, 1:arr_full.shape[1]] = ori[:, 1:ori.shape[1]]
        arr_full = arr_full[arr_full[:, 0].argsort()]

        return (arr_full)

    def add_time_dummy(self, df: DataFrame, variables: dict, _time: str):

        unique_time = sorted(df[_time].unique())
        col_timedumm = []

        prefix = _time + '_'
        for num in unique_time:
            name = prefix + str(num)
            df[name] = np.where(df[_time] == num, 1, 0)
            # new_var = regular_variable(name, 0)
            # variables['dep_indep'].append(new_var)
            # variables['iv'].append(new_var)
            col_timedumm.append(name)

        return col_timedumm

    def generate_D_matrix(self, height, T, level):
        # matrix used in Forward Orthogonal Deviation
        temp=np.zeros((T,T), dtype='float64')
        D = np.zeros((height, T), dtype='float64')

        for i in range(T):
            for j in range(i, T):
                if i == j:
                    temp[i, j] = math.sqrt((T - i - 1) / (T - i))
                else:
                    temp[i, j] = (-1) * math.sqrt(1 / ((T - i) * (T - i - 1)))

        if level:
            last=T-2


        else:
            last=T-3
        start = last + 1 - height

        D=temp[start:(last+1),:]

        return (D)
