from dataclasses import dataclass

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
    # int num_vars
    # int num_gmm_instr
