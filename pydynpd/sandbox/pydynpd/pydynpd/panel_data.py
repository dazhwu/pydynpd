import math

import numpy as np
from pandas import DataFrame


class panel_data():
    def __init__(self, df: DataFrame, identifiers):

        self._individual = identifiers[0]
        self._time = identifiers[1]

        self.N, self.T, self.ids = self.xtset(df, self._individual, self._time)
      
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
            print('Warning: system and difference GMMs do not work well on long (T>=N) panel data')
        return (N, T, ['_NT'])
    def export_data(self, df, temp_part1_list, temp_iv_list, DGMM_list, LGMM_list, timedumm):
        temp_list =  temp_part1_list.names + temp_iv_list.names + DGMM_list.names + LGMM_list.names    

        cols = list(dict.fromkeys(temp_list))

        if timedumm:
            self.col_timedumm = self.add_time_dummy(df, self._time)
        else:
            self.col_timedumm = []

        self.cols = self.ids + cols + self.col_timedumm  # make sure ids is the first column

        self.data = self.make_balanced(df[self.cols].to_numpy(), self.N, self.T)
       


    def make_balanced(self, ori, n_individual, n_time):
        arr_full = np.empty((n_individual * n_time, ori.shape[1]), dtype='float64')

        arr_full[:] = np.NaN
        
        arr_full[:, 0] = range(0, n_individual * n_time)
  

        mask = np.in1d(arr_full[:, 0], ori[:, 0])

        arr_full[mask, 1:arr_full.shape[1]] = ori[:, 1:ori.shape[1]]
        arr_full = arr_full[arr_full[:, 0].argsort()]

        return (arr_full)

    def add_time_dummy(self, df: DataFrame,  _time: str):

        unique_time = sorted(df[_time].unique())
        col_timedumm = []

        prefix = _time + '_'
        for num in unique_time:
            name = prefix + str(num)
            df[name] = np.where(df[_time] == num, 1, 0)
        
            col_timedumm.append(name)

        return col_timedumm

    