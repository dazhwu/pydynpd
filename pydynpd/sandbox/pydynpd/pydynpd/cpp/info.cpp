#include "info.h"

model_options::model_options()
{
	steps = 2;
	level = true;
	beginner = false;
	timedumm = false;
	collapse = false;
	mmsc = "bic";
	transformation = "fd";
}

basic_info::basic_info(int _N, int _T, int _num_obs, int _num_instr, int _num_indep, int _diff_width, int _max_obs, int _min_obs, double _avg_obs, vector<regular_variable> dep_indep, bool level)
{
	N = _N;
	T = _T;
	num_obs = _num_obs;
	num_instr = _num_instr;
	num_indep = _num_indep;
	diff_width = _diff_width;
	max_obs = _max_obs;
	min_obs = _min_obs;
	avg_obs = _avg_obs;
	
	for (std::size_t i = 0, max = dep_indep.size(); i < max; ++i)
	{
		if (i == 0)
			dep = dep_indep[i].name;
		else
		{
			int lag = dep_indep[i].lag;
			if (lag == 0)
				indep.push_back(dep_indep[i].name);
			else
				indep.push_back("L" + std::to_string(lag) + "." + dep_indep[i].name);
		}
	}

	if (level)
		indep.push_back("_con");
}

Regression::Regression(RowMatrixXd &reg_tab, RowMatrixXd &_vcov, RowMatrixXd &w, Ref<RowMatrixXd> _z, Ref<RowMatrixXd> _Cy, Ref<RowMatrixXd> _Cx)
{
	regression_table = reg_tab;
	vcov = _vcov;
	weighting = w;
	Z = _z;
	Cy = _Cy;
	Cx = _Cx;
}

df_info::df_info(int _N, int _T, int _f_d_i, int _f_lev_i, int _l_d_i, int _l_lev_i, int _max_lag)
{
	N = _N;
	T = _T;
	first_diff_index = _f_d_i;
	first_level_index = _f_lev_i;
	last_diff_index = _l_d_i;
	last_level_index = _l_lev_i;
	max_lag = _max_lag;
}

z_info::z_info()
{
	diff_height = 0;
	diff_width = 0;
	level_width = 0;
	level_height = 0;
	z_height = 0;
	z_width = 0;
	num_Dgmm_instr = 0;
	num_Lgmm_instr = 0;
	num_instr = 0;
}
z_info::z_info(int _diff_height, int _diff_width, int _level_width,
			   int _level_height, int _z_height, int _z_width,
			   int _num_Dgmm_instr, int _num_Lgmm_instr)
{
	diff_height = _diff_height;
	diff_width = _diff_width;
	level_width = _level_width;
	level_height = _level_height;
	z_height = _z_height;
	z_width = _z_width;
	num_Dgmm_instr = _num_Dgmm_instr;
	num_Lgmm_instr = _num_Lgmm_instr;
	num_instr = _z_height;
}
