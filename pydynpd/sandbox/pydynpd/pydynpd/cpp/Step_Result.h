//
// Created by Tiger on 5/3/2022.
//

#ifndef UNTITLED_STEP_RESULT_H
#define UNTITLED_STEP_RESULT_H

#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <tuple>

#include "Common_Functions.h"

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>
#include <eigen3/Eigen/QR>

using std::vector;
using std::string;
using std::tuple;


using Eigen::MatrixXd;
using Eigen::Map;
using Eigen::Ref;


using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

class Step_Result {
public:
    RowMatrixXd residual, _residual_t, XZ_W, W, W_inv, W_next, _M_XZ_W, zs, ZuuZ, vcov, M,
            _zs_list, beta, std_err;

    double SS;

    Step_Result(RowMatrixXd W, RowMatrixXd W_inv_, RowMatrixXd W_next_, RowMatrixXd XZ_W_,
                RowMatrixXd M_, RowMatrixXd zs_, RowMatrixXd ZuuZ_, RowMatrixXd vcov_, RowMatrixXd _M_XZ_W_,
                RowMatrixXd _zs_list_,
                RowMatrixXd beta_, RowMatrixXd std_err_, RowMatrixXd residual_, RowMatrixXd _residual_t_, double SS_);

};


#endif //UNTITLED_STEP_RESULT_H
