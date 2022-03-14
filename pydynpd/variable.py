
class regular_variable:
    def __init__(self, name, lag):
        self.name =name
        self.lag =lag


class gmm_var(regular_variable):
    def __init__(self, name, min_lag, max_lag, lag):
        super().__init__(name, lag)
        self.min_lag =min_lag
        self.max_lag =max_lag

# class indep_var(rhs_var):
#     pass
#
# class iv_var (rhs_var):
#     pass