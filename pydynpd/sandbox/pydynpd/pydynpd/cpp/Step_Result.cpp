
#include "Step_Result.h"

Step_Result::Step_Result(RowMatrixXd W_, RowMatrixXd W_inv_, RowMatrixXd W_next_, RowMatrixXd XZ_W_,
                         RowMatrixXd M_, RowMatrixXd zs_, RowMatrixXd ZuuZ_, RowMatrixXd vcov_, RowMatrixXd _M_XZ_W_,
                         RowMatrixXd _zs_list_,
                         RowMatrixXd beta_, RowMatrixXd std_err_, RowMatrixXd residual_, RowMatrixXd _residual_t_,
                         double SS_) {

    W = W_;

    W_inv = W_inv_;
    W_next = W_next_;
    XZ_W = XZ_W_;
    M = M_;
    zs = zs_;
    ZuuZ = ZuuZ_;
    vcov = vcov_;
    _M_XZ_W = _M_XZ_W_;
    _zs_list = _zs_list_;
    beta = beta_;
    std_err = std_err_;
    residual = residual_;
    _residual_t = _residual_t_;
    SS = SS_;


}

