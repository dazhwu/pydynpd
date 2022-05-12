#pragma once

#include<stdio.h>
#include<stdlib.h>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <tuple>
#include <cstring>
#include<fstream>


//#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
//#define EIGEN_USE_BLAS
//#define EIGEN_USE_LAPACKE


#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>
#include <eigen3/Eigen/QR>

using std::vector;
using std::string;
using std::tuple;


using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::Map;
using Eigen::Ref;

using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

