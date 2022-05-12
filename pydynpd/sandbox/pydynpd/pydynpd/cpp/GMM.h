//
// Created by Tiger on 5/3/2022.
//

#ifndef UNTITLED_GMM_H
#define UNTITLED_GMM_H

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include "dynpd.h"
#include "info.h"

#include "Step_Result.h"
#include "Common_Functions.h"
#include "Specification_Test.h"


//std::tuple<vector<Step_Result> , Hansen_test_info , vector<AR_test_info> >
std::tuple <struct Regression, Hansen_test_info, vector<AR_test_info>, struct basic_info> regular_process(Eigen::Ref <RowMatrixXd>,
                                                                                         Eigen::Ref <RowMatrixXd>,
                                                                                         Eigen::Ref <RowMatrixXd>,
                                                                                         struct basic_info,
				struct model_options);

RowMatrixXd regression_table(RowMatrixXd, RowMatrixXd);

void GMM(int, int, int, int, Ref <RowMatrixXd>,
         Ref <MatrixXd>,
         Ref <RowMatrixXd>,
         Ref <RowMatrixXd>, Ref <RowMatrixXd>,
         Ref <RowMatrixXd>, Ref <RowMatrixXd>,
         Ref <RowMatrixXd>, int step, string transformation, bool level);

std::tuple <RowMatrixXd, RowMatrixXd> calculate_basic(Ref <RowMatrixXd> z_list,
                                                      Ref <RowMatrixXd> Cx_list,
                                                      Ref <RowMatrixXd> Cy_list, int N);

RowMatrixXd calculate_W(Ref <RowMatrixXd> H,
                        Ref <RowMatrixXd> z_list,
                        Ref <MatrixXd> _z_t_list, int step,
                        int N);

RowMatrixXd calculate_residual(Ref <RowMatrixXd> y_list,
                               Ref <RowMatrixXd> x_list,
                               Ref <RowMatrixXd> beta);

RowMatrixXd vcov_step_1(Ref <RowMatrixXd> _M_XZ_W, Ref <RowMatrixXd> W2, int N);

RowMatrixXd vcov(Ref <RowMatrixXd> z_list,
                 Ref <RowMatrixXd> Cx, Ref <RowMatrixXd> M2, Ref <RowMatrixXd> _M2_XZ_W2, Ref <RowMatrixXd> _W2_inv,
                 Ref <RowMatrixXd> zs2,
                 int step, int N);


Hansen_test_info perform_hansen_test(int step, int num_instru, int num_indep, int N);

vector <AR_test_info> perform_AR_test(Ref <RowMatrixXd> Cx, Ref <RowMatrixXd> Cy,
                                      Ref <RowMatrixXd> z_list, Ref <RowMatrixXd> _zs_list,
                                      int step, int diff_width, int N, string transformation, bool level);

/*
std::tuple<Hansen_test_info, vector<AR_test_info>>  perform_test( Ref<RowMatrixXd> ,  Ref<RowMatrixXd> ,
                     Ref<RowMatrixXd> , Ref<RowMatrixXd> ,
                  int num_instru, int num_indep, int step, int didd_width,	int N, string transformation, bool level);
*/
RowMatrixXd get_H1(Ref <RowMatrixXd> z_list, int diff_width, int T,
                   string transformation, bool level);

RowMatrixXd get_H1_fod(int width, int diff_width, int T, bool level);

RowMatrixXd generate_D_matrix(int height, int T);

RowMatrixXd Windmeijer(Ref <RowMatrixXd> M2, Ref <RowMatrixXd> _M2_XZ_W2,
                       Ref <RowMatrixXd> W2_inv, Ref <RowMatrixXd> zs2,
                       Ref <RowMatrixXd> vcov_step1, Ref <RowMatrixXd> Cx_list,
                       Ref <RowMatrixXd> z_list, Ref <RowMatrixXd> residual1_t, int N);


#endif //UNTITLED_GMM_H