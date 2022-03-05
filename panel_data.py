import pandas as pd
import numpy as np
from pandas import DataFrame

class pandel_data(object):


#https://stackoverflow.com/questions/29352511/numpy-sort-ndarray-on-multiple-columns

# https://itecnote.com/tecnote/python-efficiently-applying-a-function-to-a-grouped-pandas-dataframe-in-parallel/

    def __init__(self, df: DataFrame, identifiers):
        pass
        # turn df to numpy array


    def make_balanced(df: object, N, T) -> object:

        N=5
        T=2
        left=np.array([[0, 0, 0, 2],[1, 0, 1, 12]])
        arr_full =np.empty((N*T, 4), dtype='float64')
        arr_full[:] = np.NaN
        arr_full[:,0]=range(0, N*T)

        mask = np.in1d(arr_full[:, 0],left[:, 0])
        arr_full[mask, 1:5]=left[:,1:4]

