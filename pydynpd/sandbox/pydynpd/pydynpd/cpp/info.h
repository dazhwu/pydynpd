#pragma once



#include "dynpd.h"
#include "List_Variables.h"
#include "variable.h"



struct model_options {
  int steps;             // = 2;
  bool level;            // = true;
  bool beginner;         // = false;
  bool timedumm;         // = false;
  bool collapse;         // = false;
  string mmsc;           // = "bic";
  string transformation; // = "fd";
  model_options();
};


struct df_info {
	int N;
	int T;
	//vector<string> ids;
	int first_diff_index;
	int last_diff_index;
	int first_level_index;
	int last_level_index;
	int max_lag;
	df_info(int _N, int _T, int _f_d_i, int _f_lev_i, int _l_d_i, int _l_lev_i, int _max_lag);
	// last_fod_index;
	//# first_fod_index; int
};

struct z_info {
  int diff_width;
  int diff_height;
  int level_width;
  int level_height;
  int z_width;
  int z_height;
  int num_Dgmm_instr;
  int num_Lgmm_instr;
  int num_instr;
  z_info();
  z_info(int, int, int, int, int, int, int, int);
  //# int num_vars
  //# int num_gmm_instr
};



struct Regression {
   RowMatrixXd regression_table;
   RowMatrixXd vcov;
   RowMatrixXd weighting;
   RowMatrixXd Z;
   RowMatrixXd Cy;
   RowMatrixXd Cx;
   
   
  Regression(RowMatrixXd &reg_tab, RowMatrixXd &_vcov, RowMatrixXd &w, Ref<RowMatrixXd> _z, Ref<RowMatrixXd> _Cy, Ref<RowMatrixXd> _Cx);
};


struct basic_info {
    string dep;
    vector<string> indep;
    int N;
    int T;
    int num_obs;
    int num_instr;
    int num_indep;
    int diff_width;
    int max_obs;
    int min_obs;
    double avg_obs;
    double SS;
    basic_info(int, int, int, int,  int,int,int,int,double, vector<regular_variable>, bool);
};


