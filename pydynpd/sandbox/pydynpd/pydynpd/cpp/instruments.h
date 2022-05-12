#pragma once
#include "dynpd.h"
#include "info.h"
#include "List_Variables.h"
#include "variable.h"


class instruments {
  struct z_info z_information;
  RowMatrixXd z_table;
};

std::tuple<RowMatrixXd, z_info>
get_z_table(int N, int T, struct df_info info, Ref<RowMatrixXd> Dgmm_dat,
            Ref<RowMatrixXd> Lgmm_dat, Ref<RowMatrixXd> iv_dat,
            Ref<RowMatrixXd> Div_dat, vector<gmm_var> Dgmm_vars,
            vector<gmm_var> Lgmm_vars, vector<regular_variable>iv_vars, bool level,
            string transformation, bool collapse);



void build_z_level(int N, vector<gmm_var> Dgmm_vars, vector<gmm_var> Lgmm_vars,
                   vector<regular_variable> iv_vars, struct df_info info,
              struct z_info z_information, Ref<RowMatrixXd> Dgmm_dat,
              Ref<RowMatrixXd> Lgmm_dat, Ref<RowMatrixXd> Div_dat,
              Ref<RowMatrixXd> iv_dat, bool level, string transformation,
              bool collapse) ;

int prepare_z_gmm_level(vector<gmm_var> Lgmm_vars, struct df_info info, int level_width,
                        bool collapse ) ;

struct z_info calculate_z_dimension(vector<gmm_var> Dgmm_vars,
                                    vector<gmm_var> Lgmm_vars,
                                    vector<regular_variable> iv_vars, struct df_info info,
                                    bool level, string transformation,
                                    bool collapse) ;

void build_z_diff(int N, vector<gmm_var> Dgmm_vars, vector<regular_variable> iv_vars,
                  struct df_info info, struct z_info z_information,
                  Ref<RowMatrixXd> Dgmm_dat, Ref<RowMatrixXd> Div_dat,
                  bool level, string transformation, bool collapase ) ;

void prepare_Z_iv_diff(vector<regular_variable> iv_vars, int width, struct df_info info) ;

int prepare_Z_gmm_diff(vector<gmm_var> Dgmm_vars, struct df_info info,
                       bool level, string transformation,
                       bool collapse = false);