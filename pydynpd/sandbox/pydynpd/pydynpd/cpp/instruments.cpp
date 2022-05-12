#include "instruments.h"

#define EIGEN_INITIALIZE_MATRICES_BY_NAN
RowMatrixXd z_table;
Eigen::MatrixXi gmm_diff_info, iv_diff_info, gmm_level_info;

// int level_width, level_height, diff_width, diff_height, z_width, z_height;
// int Lgmm_iv_height, Lgmm_iv_width;


std::tuple <RowMatrixXd, z_info>
get_z_table(int N, int T, struct df_info info, Ref <RowMatrixXd> Dgmm_dat,
            Ref <RowMatrixXd> Lgmm_dat, Ref <RowMatrixXd> iv_dat,
            Ref <RowMatrixXd> Div_dat, vector <gmm_var> Dgmm_vars,
            vector <gmm_var> Lgmm_vars, vector <regular_variable> iv_vars, bool level,
            string transformation, bool collapse) {
	struct z_info z_information = calculate_z_dimension(
			Dgmm_vars, Lgmm_vars, iv_vars, info, level, transformation, collapse);
	/*
	level_width = z_information.level_width;
	level_height = z_information.level_height;
	diff_width = z_information.diff_width;
	diff_height = z_information.diff_height;
	z_width = z_information.z_width;
	z_height = z_information.z_height;
	*/
	if (level)
		build_z_level(N, Dgmm_vars, Lgmm_vars, iv_vars, info, z_information,
		              Dgmm_dat, Lgmm_dat, Div_dat, iv_dat, level, transformation,
		              collapse);
	else
		build_z_diff(N, Dgmm_vars, iv_vars, info, z_information, Dgmm_dat, Div_dat,
		             level, transformation, collapse);

    z_table = z_table.unaryExpr([](double v) { return std::isfinite(v)? v : 0.0; });
	return std::make_tuple(z_table, z_information);
}

void build_z_level(int N, vector <gmm_var> Dgmm_vars, vector <gmm_var> Lgmm_vars,
                   vector <regular_variable> iv_vars, struct df_info info,
                   struct z_info z_information, Ref <RowMatrixXd> Dgmm_dat,
                   Ref <RowMatrixXd> Lgmm_dat, Ref <RowMatrixXd> Div_dat,
                   Ref <RowMatrixXd> iv_dat, bool level, string transformation,
                   bool collapse) {
	int level_info_width = gmm_level_info.cols();

	build_z_diff(N, Dgmm_vars, iv_vars, info, z_information, Dgmm_dat, Div_dat,
	             level, transformation, collapse);
				 

	int start_row = z_information.diff_height; //# z_table[0].shape[0] - height
	int start_col = z_information.diff_width;  //# z_table[0].shape[1] - width

	int Lgmm_iv_height = Lgmm_dat.rows() / N;
	int Lgmm_width = Lgmm_dat.cols();
	int iv_width=iv_dat.cols();
	//std::cout << start_col << " " << z_information.z_width<< std::endl;
	//std::cout << start_row << " " << z_information.z_height << std::endl;
	//Eigen::MatrixXi temp_col=Eigen::MatrixXi::Constant(1,1,  z_information.z_width-start_col);
	#pragma omp parallel for
	for (int i = 0; i < N; ++i) {
		Ref <RowMatrixXd> z =
				z_table.block(i * z_information.z_height, 0, z_information.z_height,
				              z_information.z_width);
		//# z[start_row - 1, start_col:(start_col + self.level_width)] = 1
		for (int k = start_col; k < z_information.z_width; ++k)
			z(z_information.z_height - 1, k) = 1; // part of last row


		Ref <RowMatrixXd> array_Lgmm = Lgmm_dat.block(i * Lgmm_iv_height, 0, Lgmm_iv_height, Lgmm_width);
		Ref <RowMatrixXd> array_iv = iv_dat.block(i * Lgmm_iv_height, 0, Lgmm_iv_height, iv_width);
		for (std::size_t var_id = 0, max = Lgmm_vars.size(); var_id != max; ++var_id) {
			//int lag = Lgmm_vars.min_lags[var_id];

			Ref <MatrixXi> t_info =
					gmm_level_info.block(3 * var_id, 0, 3, level_info_width);
			for (int j = 0; j < z_information.level_width; ++j) {
				int the_index = t_info(0, j);
				if (the_index >= 0) {
					int the_row = start_row + t_info(2, j);
					// the_index = t_info[0, j];
					z(the_row, start_col + j) = array_Lgmm(the_index, var_id);
					
				}
			}
		}
		
		int start_pos = z_information.num_Dgmm_instr;
		int len_iv_vars = iv_vars.size();
		for (int var_id = 0; var_id < len_iv_vars; ++var_id) {
			for (int j = 0; j < z_information.level_width; ++j) {
				z(start_pos + var_id, j + z_information.diff_width) =
						array_iv(info.first_level_index + j, var_id);
			}
		//std::cout << the_row << " " << start_col+j << array_Lgmm(the_index, var_id) << std::endl;
		}
	}

	//   z[np.isnan(z)] = 0
}

int prepare_z_gmm_level(vector <gmm_var> Lgmm_vars, struct df_info info, int level_width,
                        bool collapse) {
	int num_gmm = Lgmm_vars.size();
	gmm_level_info = MatrixXi::Zero(num_gmm * 3, level_width);
	// t_info = np.empty((num_gmm * 3, level_width),                  dtype =
	// 'int32')  # row 0: gmm_index, 2 : which row, 1 : empty right now
	
	int start_row = 0;
	for (int var_id = 0; var_id < num_gmm; ++var_id) {
		// var = Lgmm_vars[var_id]
		for (int i = 0; i < level_width; ++i) {
			int the_index = info.first_level_index + i;
			int gmm_index = the_index - Lgmm_vars[var_id].min_lag;
			gmm_level_info(var_id * 3 + 2, i) = start_row;
			if (gmm_index >= 1){
				gmm_level_info(var_id * 3 + 0, i) = gmm_index;
				if (collapse) {
					if (i == (level_width - 1))
						start_row += 1;
				}else
					start_row += 1;
			}else
				gmm_level_info(var_id * 3 + 0, i) = -9999; // need to change
		}
	}

	
	// num_Lgmm_instr = start_row; //#number of gmm instruments in diff eq
	return start_row;
}

struct z_info calculate_z_dimension(vector <gmm_var> Dgmm_vars,
                                    vector <gmm_var> Lgmm_vars,
                                    vector <regular_variable> iv_vars, struct df_info info,
                                    bool level, string transformation,
                                    bool collapse) {
	// Lgmm_vars = variables['Lgmm']

	int diff_width = info.last_diff_index - info.first_diff_index + 1;
	int level_width = info.last_level_index - info.first_level_index + 1;

	int level_height = 0;
	int num_Lgmm_instr = 0;

	// gmm_level_info = None
	/*void prepare_z_gmm_level(vector<gmm_var> Lgmm_vars, variables
						   : dict, struct  df_info info,
							 bool collapse = False) */
	if (level) {
		int num_Lgmm_instr = prepare_z_gmm_level(Lgmm_vars, info, level_width, collapse);
		level_height = num_Lgmm_instr + 1;
	}
	int num_Dgmm_instr =
			prepare_Z_gmm_diff(Dgmm_vars, info, level, transformation, collapse);
	prepare_Z_iv_diff(iv_vars, diff_width, info);
	int diff_height = (num_Dgmm_instr + iv_diff_info.rows());
	int z_height, z_width;
	if (level) {
		z_height = diff_height + level_height;
		z_width = diff_width + level_width;
	} else {
		z_height = diff_height;
		z_width = diff_width;
	}

	

	z_info z_information =
			z_info(diff_height, diff_width, level_width, level_height, z_height,
			       z_width, num_Dgmm_instr, num_Lgmm_instr);
	return z_information; //   , gmm_diff_info, iv_diff_info, gmm_level_info)
}

void build_z_diff(int N, vector <gmm_var> Dgmm_vars, vector <regular_variable> iv_vars,
                  struct df_info info, struct z_info z_information,
                  Ref <RowMatrixXd> Dgmm_dat, Ref <RowMatrixXd> Div_dat,
                  bool level, string transformation, bool collapase) {
	int gmm_Div_height = Dgmm_dat.rows() / N;
	int gmm_width = Dgmm_dat.cols();
	int div_width =Div_dat.cols();

	int diff_width = z_information.diff_width;
	int z_height = z_information.z_height;
	int z_width = z_information.z_width;
	int num_Dgmm_instr = z_information.num_Dgmm_instr;
	
	z_table = RowMatrixXd::Zero(z_height * N, z_width);
	int start_col = 0;
	if ((level) && (transformation == "fod"))
		start_col = 1;
	#pragma omp parallel for
	for (int i = 0; i < N; ++i) {
		Ref <RowMatrixXd> z = z_table.block(i * z_height, start_col, z_height, diff_width - start_col);

		// z_width = z.shape[1]
		Ref <RowMatrixXd> array_gmm =
				Dgmm_dat.block(i * gmm_Div_height, 0, gmm_Div_height, gmm_width);
				
				
		Ref <RowMatrixXd> array_fd_iv =
				Div_dat.block(i * gmm_Div_height, 0, gmm_Div_height, div_width);
				
		for (std::size_t var_id = 0, max = Dgmm_vars.size(); var_id != max; ++var_id) {
			for (int j = 0; j < z.cols(); ++j) {
				int row_pos = gmm_diff_info(var_id * 3 + 2, j);

				int start = gmm_diff_info(var_id * 3 + 0, j);
				int end = gmm_diff_info(var_id * 3 + 1, j);
			
				for (int k = 0; k < (end - start + 1); ++k)
					z(row_pos + k, j) = array_gmm(end - k, var_id);
			
			}
			
		}
		
		for (std::size_t var_id = 0, max = iv_vars.size(); var_id != max; ++var_id) {
			int row_pos = num_Dgmm_instr + var_id;
			for (int j = 0; j < z.cols(); ++j) {
				int index_to_take = iv_diff_info(var_id, j);
				z(row_pos, j) = array_fd_iv(index_to_take, var_id);
			}
		}
		
	}

	//saveData("z_tab.csv", z_table);
}

void prepare_Z_iv_diff(vector <regular_variable> iv_vars, int width, struct df_info info) {

	// iv_vars = variables['iv'] #need to be placed at the beginning

	int num_iv = iv_vars.size();
	iv_diff_info = Eigen::MatrixXi::Zero(num_iv, width);
	int j;
	for (j = 0; j < width; ++j) {
		Eigen::VectorXi tempCol = Eigen::VectorXi::Constant(num_iv, 1, info.first_diff_index + j);
		iv_diff_info.col(j) = tempCol;
	}

	
	/*
	for (int var_id = 0; var_id < num_iv; ++var_id)
	  for (int j = 0; j < width; ++j)
		iv_diff_info(var_id, j) = info.first_diff_index + j;*/
}

int prepare_Z_gmm_diff(vector <gmm_var> Dgmm_vars, struct df_info info,
                       bool level, string transformation,
                       bool collapse) {
	int num_gmm = Dgmm_vars.size();
	int first_index;
	if ((level) && (transformation == "fod"))
		first_index = info.first_diff_index + 1;
	else
		first_index = info.first_diff_index;
	int last_index = info.last_diff_index;
	int width = last_index - first_index + 1;
	gmm_diff_info = MatrixXi::Zero(num_gmm * 3, width);
	int start_row = 0;
	for (int var_id = 0; var_id < num_gmm; ++var_id) {
		//var = gmm_vars[var_id]

		int first_tend = first_index - Dgmm_vars[var_id].min_lag;
		//int last_tend = last_index - Dgmm_vars.min_lags[var_id];
		for (int j = 0; j < width; ++j) {
			int tend = first_tend + j;
			gmm_diff_info(var_id * 3 + 1, j) = tend;
			int tstart =
					tend + Dgmm_vars[var_id].min_lag - Dgmm_vars[var_id].max_lag;
			if (tstart < 0)
				tstart = 0;
			gmm_diff_info(var_id * 3, j) = tstart;
			gmm_diff_info(var_id * 3 + 2, j) = start_row;
			int num = tend - tstart + 1;
			//#num_instru += num
			if (collapse) {
				if (j == (width - 1))
					start_row += num;

			} else
				start_row += num;
		}
	}
	
	// num_gmm_instr = start_row // # number of gmm instruments in diff eq
	return start_row;
}

