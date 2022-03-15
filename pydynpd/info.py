from dataclasses import dataclass

import numpy as np


@dataclass
class df_info:
    N: int
    T: int
    _individual: str
    _time: str
    ids: list
    first_index: int
    last_index: int
    max_lag: int



@dataclass
class z_info:
    diff_width: int
    diff_height: int
    level_width: int
    level_height: int
    num_gmm_instr: int
    num_instr: int
    # int num_vars
    # int num_gmm_instr

@dataclass
class hansen_test_info:
    test_value: float
    df: int
    p_value: float
    critical_value: float

@dataclass
class AR_test_info:
    lag: int

@dataclass
class options_info:
    twosteps: bool=True
    level: bool=True
    timedumm: bool=False
    collapse: bool=False


# @dataclass
# class regression_info:
#     twosteps: bool
#     robust: bool
#     level: bool
#     num_obs: int
#     num_instruments: int
#     N:int
#     T: int
#     dep: str
#     indep: list
#     beta1: list
#     beta2: list
#     residual1: list
#     residual1_t: list
#     SS: float
#     XZ_W1: np.ndarray
#     XZ_W2: np.ndarray
#     W1: np.ndarray
#     W2: np.ndarray
#     M1: np.ndarray
#     M2: np.ndarray
#
#     def __init__(self, twosteps, robust, level, num_obs,  num_instruments, N, T,
#                  beta1, residual1, residual1_t, SS, XZ_W1,  XZ, M1, M2, W1, W2):
#         self.twosteps=twosteps
#         self.robust=robust
#         self.level=level
#         self.num_obs=num_obs
#         self.num_instruments=num_instruments
#         self.N=N
#         self.T=T
#         self.beta1=beta1
#         self.residual1=residual1
#         self.residual1_t=residual1_t
#         self.SS=SS
#         self.XZ_W1=XZ_W1
#         self.XZ=XZ
#         self.M1=M1
#         self.M2 = M2
#         self.W1=W1
#         self.W2=W2
#
