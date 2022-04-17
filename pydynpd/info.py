from dataclasses import dataclass

# from pandas import DataFrame
import numpy as np


@dataclass
class df_info:
    N: int
    T: int
    ids: list
    first_diff_index: int
    last_diff_index: int
    first_level_index: int
    last_level_index: int
    max_lag: int
    # last_fod_index: int
    # first_fod_index: int


@dataclass
class z_info:
    diff_width: int
    diff_height: int
    level_width: int
    level_height: int
    width: int
    height: int
    num_Dgmm_instr: int
    num_Lgmm_instr: int
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
    AR: float
    P_value: float


@dataclass
class options_info:
    steps: int = 2
    level: bool = True
    beginner: bool = False
    timedumm: bool = False
    collapse: bool = False
    mmsc: str = 'bic'
    transformation: str = 'fd'


@dataclass
class sumproduct_task:
    array_list: list
    division_list: list


@dataclass
class beginner_models:
    model: str


@dataclass
class step_result:
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
    std_err: list

    def __init__(self, W):
        self.W = W
        self.W_inv = np.linalg.pinv(W)
