#pragma once

#include "dynpd.h"
#include "Common_Functions.h"
#include <boost/math/distributions/chi_squared.hpp>

class AR_test_info {
public:
    int lag;
    double AR;
    double P_value;

    AR_test_info(int _lag, double _AR, double _P);

};

class Hansen_test_info {
public:
    double test_value;
    int df;
    double P_value;
    double critical_value;

    Hansen_test_info(double, int, double, double);

};


Hansen_test_info hansen_overid(const Ref<const RowMatrixXd> &, const Ref<const RowMatrixXd> &, int, int, int);

vector<MatrixXd> AR_get_diff_XR(int N, Eigen::Ref <RowMatrixXd> beta, Eigen::Ref <RowMatrixXd> ori_residual,
                                      Eigen::Ref <RowMatrixXd> ori_x, int diff_width,
                                      Eigen::Ref <RowMatrixXd> Diff_y_table, Eigen::Ref <RowMatrixXd> Diff_x_table,
                                      int r0_height, string transformation, bool level);

vector<AR_test_info> AR_test(int N, int m, Eigen::Ref <RowMatrixXd> z_list, Eigen::Ref <RowMatrixXd> zs_list,
                                   Eigen::Ref <RowMatrixXd> ori_residual,
                                   Eigen::Ref <RowMatrixXd> M_XZ_W, Eigen::Ref <RowMatrixXd> vcov,
                                   Eigen::Ref <RowMatrixXd> ori_x, Eigen::Ref <RowMatrixXd> beta,
                                   int diff_width, Eigen::Ref <RowMatrixXd> Diff_y_table,
                                   Eigen::Ref <RowMatrixXd> Diff_x_table, string transformation, bool level);

