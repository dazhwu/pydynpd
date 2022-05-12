#pragma once
#include <unordered_map>
#include <math.h>
#include "dynpd.h"
#include "info.h"
#include "variable.h"
#include "instruments.h"
#include "Common_Functions.h"
#include "Command.h"
#include "GMM.h"




struct df_info get_info(model_options options);

std::tuple  <struct Regression, Hansen_test_info, vector<AR_test_info>, struct basic_info>
prepare_data(Ref <RowMatrixXd> pdata,  List_Variables list_dep_indep, List_Variables list_iv,
                  List_Variables list_gmm,                   model_options options, vector<string>, vector<int>,  vector<string>, vector<string>);

void get_gmm_tables(Ref<RowMatrixXd> pdata, Ref<RowMatrixXd> fd_data,vector <gmm_var> Dgmm_vars, vector <regular_variable> iv_vars,
                    vector <gmm_var> Lgmm_vars, vector<string>,bool level) ;
template<class ty>
RowMatrixXd gen_table(Ref <RowMatrixXd> ori_data, vector <ty> variable_list,
                      vector <string> df_cols) ;

void get_xy_table_dict(Ref <RowMatrixXd> pdata, vector <regular_variable> dep_indep, vector<string>,
string transformation);

void get_final_xy_tables(df_info info, bool level, string transformation) ;

std::tuple <RowMatrixXd, RowMatrixXd> get_final_xy_systemGMM(
		Ref <RowMatrixXd> Dx, Ref <RowMatrixXd> Dy, Ref <RowMatrixXd> x, Ref <RowMatrixXd> y,
		vector<int> cut, vector<int> Dcut, int cut_height, int Dcut_height, string transformation) ;



std::tuple <RowMatrixXd, RowMatrixXd>
get_final_xy_diffGMM(Ref <RowMatrixXd> Dy, Ref <RowMatrixXd> Dx, vector<int> Dcut, int Dcut_height);

void prepare_reg_fod(Ref <RowMatrixXd> Diff_x, Ref <RowMatrixXd> Diff_y);

std::tuple<int, int, int, double>
prepare_reg(Ref <RowMatrixXd> z_list, Ref <RowMatrixXd> Cx, Ref <RowMatrixXd> Cy,z_info,
             string transformation, bool level);