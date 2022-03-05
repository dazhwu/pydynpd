from dataclasses import dataclass

@dataclass
class dep_var:
    name: str

@dataclass
class indep_var:
    name: str
    lag_start: int
    lag_end: int

@dataclass
class gmm_var:
    name: str
    min_lag: int
    max_lag: int
    diff_level: int

@dataclass
class iv_var:
    name: str
    lag_start: int
    lag_end: int
    diff_level: int
