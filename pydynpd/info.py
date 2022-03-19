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
    width:int
    height: int
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
    steps: int=2
    level: bool=True
    timedumm: bool=False
    collapse: bool=False


@dataclass
class sumproduct_task:
    array_list: list
    division_list: list


@dataclass
class regression_info:
    # twosteps: bool
    # robust: bool
    # level: bool
    # num_obs: int
    # num_instruments: int
    # N:int
    # T: int
    # # dep: str
    # indep: list

    # beta2: list
    H: np.ndarray
    residual: np.ndarray
    residual_t: np.ndarray
    SS: float
    XZ: np.ndarray
    Zy: np.ndarray
    XZ_W: np.ndarray
    W_inv: np.ndarray
    M_XZ_W: np.ndarray
    ZuuZ: np.ndarray
    W: np.ndarray
    M: np.ndarray
    beta: list
    vcov: np.ndarray
    def __init__(self, W):
        self.W=W
        self.W_inv=np.linalg.pinv(W)
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
