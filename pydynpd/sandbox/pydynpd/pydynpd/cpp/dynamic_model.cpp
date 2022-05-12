#include "dynamic_model.h"

int N, T;

RowMatrixXd Dgmm_table, iv_table, Delta_iv_table, Lgmm_table;
RowMatrixXd ori_x_table, ori_y_table, ori_Diff_x_table, ori_Diff_y_table;
vector<gmm_var> Dgmm_vars, Lgmm_vars;
vector<regular_variable> dep_indep, iv_vars;
std::unordered_map<string, RowMatrixXd> xy_tables, final_xy_tables;


struct df_info get_info(model_options options)
{

	int max_lag = 0;
	int max_Dgmm_minlag = 0;
	int max_Lgmm_minlag = 0;

	for (auto var : dep_indep)
	{
		if (var.lag > max_lag)
			max_lag = var.lag;
	}

	for (auto var : iv_vars)
	{
		if (var.lag > max_lag)
			max_lag = var.lag;
	}

	for (auto var : Dgmm_vars)
	{
		if (var.min_lag > max_Dgmm_minlag)
			max_Dgmm_minlag = var.min_lag;
	}

	for (auto var : Lgmm_vars)
	{
		if (var.min_lag > max_Lgmm_minlag)
			max_Lgmm_minlag = var.min_lag;
	}

	int last_level_index = T - 1;
	int last_diff_index = T - 1;

	int first_level_index, first_diff_index;
	if (max_lag > max_Lgmm_minlag)
		first_level_index = max_lag;
	else
		first_level_index = max_Lgmm_minlag;

	if ((options.level) && (options.transformation == "fod"))
		first_diff_index = first_level_index;
	else
		first_diff_index = first_level_index + 1; //# max(max_lag + 1, max_Dgmm_minlag)

	if (first_diff_index + 2 > last_diff_index) // # to do: change 3 to something rated to AR(p)
		throw std::invalid_argument("Not enough periods to run the model");

	struct df_info tbr = df_info(N, T, first_diff_index, first_level_index, last_diff_index, last_level_index, max_lag);

	return tbr;
}

std::tuple<struct Regression, Hansen_test_info, vector<AR_test_info>, struct basic_info>
prepare_data(Ref<RowMatrixXd> pdata, List_Variables list_dep_indep, List_Variables list_iv,
			 List_Variables list_gmm, model_options options, vector<string> parts, vector<int> dimensions, vector<string> df_cols, vector<string> col_timedumm)
{

	N = dimensions[0];
	T = dimensions[1];

	xy_tables.clear();
	final_xy_tables.clear();
	RowMatrixXd fd_data = get_first_diff_table(pdata, N);
	std::tie(Dgmm_vars, Lgmm_vars) = list_to_Gmm(list_gmm);
	dep_indep = list_to_dep_indep_iv(list_dep_indep); 
	iv_vars = list_to_dep_indep_iv(list_iv);

	get_gmm_tables(pdata, fd_data, Dgmm_vars, iv_vars, Lgmm_vars, df_cols, options.level);

	get_xy_table_dict(pdata, dep_indep, df_cols, options.transformation);

	struct df_info info = get_info(options);

	get_final_xy_tables(info, options.level, options.transformation);

	RowMatrixXd z_table;

	struct z_info z_information;
	std::tie(z_table, z_information) =
		get_z_table(N, T, info, Dgmm_table, Lgmm_table, iv_table, Delta_iv_table, Dgmm_vars, Lgmm_vars, iv_vars, options.level, options.transformation, options.collapse);

	int num_obs, max_obs, min_obs;
	double avg_obs;
	std::tie(num_obs, max_obs, min_obs, avg_obs) = prepare_reg(z_table, final_xy_tables["Cx"], final_xy_tables["Cy"], z_information, options.transformation,
															   options.level);
	


	int num_indep = final_xy_tables["Cx"].cols();



	struct basic_info model_info(N, T, num_obs, z_information.num_instr, num_indep, z_information.diff_width, max_obs, min_obs, avg_obs, dep_indep, options.level);

	return regular_process(z_table, final_xy_tables["Cx"], final_xy_tables["Cy"], model_info, options);
}

void get_gmm_tables(Ref<RowMatrixXd> pdata, Ref<RowMatrixXd> fd_data, vector<gmm_var> Dgmm_vars, vector<regular_variable> iv_vars,
					vector<gmm_var> Lgmm_vars, vector<string> df_cols, bool level)
{
	Dgmm_table = gen_table(pdata, Dgmm_vars, df_cols);
	iv_table = gen_table(pdata, iv_vars, df_cols);
	Delta_iv_table = gen_table(fd_data, iv_vars, df_cols);

	if (level)
		Lgmm_table = gen_table(fd_data, Lgmm_vars, df_cols);
}

template <class ty>
RowMatrixXd gen_table(Ref<RowMatrixXd> ori_data, vector<ty> variable_list,
					  vector<string> df_cols)
{
	int num_variables = variable_list.size();
	//        list_cols = self.pdata.ids.copy();

	// variable_names = vairable_list.names; //   [var.name for var in variable_list]
	vector<int> which_col;
	for (auto var : variable_list)
	{
		int the_index = getIndex(df_cols, var.name);
		which_col.push_back(the_index);
	}

	//         VectorXd A;
	// A = (A.array().isfinite()).select(A,0);
	// start_row = 0
	// end_row = T - 1

	int height = T; // end_row - start_row + 1
	RowMatrixXd tbr(height * N, num_variables);
	// int ori_width = ori_data.cols();
	for (int j = 0; j < num_variables; ++j)
	{
		int lag = variable_list[j].lag;
		if (lag == 0)
			tbr.col(j) =
				ori_data.col(which_col[j]);
		else
		{
			for (int i = 0; i < N; ++i)
			{
				// int col = 0;
				Ref<RowMatrixXd> ori_i = ori_data.block(i * T, which_col[j], T, 1);
				Ref<RowMatrixXd> tbr_i = tbr.block(i * height, j, height, 1);
				tbr_i.block(0, 0, lag, 1) = RowMatrixXd(lag, 1);
				tbr_i.block(lag, 0, height - lag, 1) =
					ori_i.block(0, 0, height - lag, 1);
			}
		}
	}
	return tbr;
}

void get_xy_table_dict(Ref<RowMatrixXd> pdata, vector<regular_variable> dep_indep, vector<string> df_cols,
					   string transformation)
{
	int num_var = dep_indep.size();
	vector<regular_variable> dep{&dep_indep[0], &dep_indep[1]};
	vector<regular_variable> indeps{&dep_indep[1], &dep_indep[num_var]};
	RowMatrixXd ori_y_table = gen_table(pdata, dep, df_cols);
	RowMatrixXd ori_x_table = gen_table(pdata, indeps, df_cols);

	/*	RowMatrixXd gen_table(Ref <RowMatrixXd> ori_data, vector <ty> variable_list,
							  vector <string> df_cols)
	*/
	xy_tables["x"] = ori_x_table;
	xy_tables["y"] = ori_y_table;
	RowMatrixXd ori_Diff_y_table = get_first_diff_table(ori_y_table, N);
	RowMatrixXd ori_Diff_x_table = get_first_diff_table(ori_x_table, N);
	if (transformation == "fd")
	{
		xy_tables["Dy"] = ori_Diff_y_table;
		xy_tables["Dx"] = ori_Diff_x_table;
	}
	else
	{
		// ori_Fod_y_table = get_fod_table(ori_y_table, N);
		// ori_Fod_x_table = get_fod_table(ori_x_table, N);
		xy_tables["Dy"] = get_fod_table(ori_y_table, N);
		xy_tables["Dx"] = get_fod_table(ori_x_table, N);
		xy_tables["Diff_y"] = ori_Diff_y_table;
		xy_tables["Diff_x"] = ori_Diff_x_table;
	}
}

void get_final_xy_tables(df_info info, bool level, string transformation)
{
	vector<int> Dcut = {info.first_diff_index, info.last_diff_index};
	vector<int> cut = {info.first_level_index, info.last_level_index};

	int Dcut_height = Dcut[1] - Dcut[0] + 1;
	int cut_height = cut[1] - cut[0] + 1;

	RowMatrixXd Cy, Cx, Diff_y, Diff_x;

	if (level)
	{
		std::tie(Cy, Cx) = get_final_xy_systemGMM(
			xy_tables["Dx"], xy_tables["Dy"], xy_tables["x"], xy_tables["y"], cut, Dcut, cut_height, Dcut_height, transformation);
	}
	else
	{
		std::tie(Cy, Cx) = get_final_xy_diffGMM(
			xy_tables["Dy"], xy_tables["Dx"], Dcut, Dcut_height);
	}
	final_xy_tables["Cy"] = Cy;
	final_xy_tables["Cx"] = Cx;

	if (transformation == "fod")
	{
		std::tie(Diff_y, Diff_x) = get_final_xy_diffGMM(
			xy_tables["Diff_y"], xy_tables["Diff_x"], Dcut, Dcut_height);
		if (level)
		{
			int height = Diff_y.rows();
			RowMatrixXd zeros = RowMatrixXd::Zero(height, 1);

			Diff_x.conservativeResize(height, Diff_x.cols() + 1);

			Diff_x.col(Diff_x.cols() - 1) = zeros;
			// Diff_x = np.hstack((Diff_x, zeros))
		}
		final_xy_tables["Diff_y"] = Diff_y;
		final_xy_tables["Diff_x"] = Diff_x;
	}
}

std::tuple<RowMatrixXd, RowMatrixXd> get_final_xy_systemGMM(
	Ref<RowMatrixXd> Dx, Ref<RowMatrixXd> Dy, Ref<RowMatrixXd> x, Ref<RowMatrixXd> y,
	vector<int> cut,
	vector<int> Dcut, int cut_height,
	int Dcut_height, string transformation)
{

	int height_total = Dcut_height + cut_height;

	int width = x.cols();
	int Dx_height = Dx.rows() / N;
	int x_height = x.rows() / N;

	RowMatrixXd Cx(height_total * N, width + 1);
	RowMatrixXd Cy(height_total * N, 1);
#pragma omp parallel for
	for (int i = 0; i < N; ++i)
	{
		Ref<RowMatrixXd> Dy_i = Dy.block(i * Dx_height, 0, Dx_height, 1);
		Ref<RowMatrixXd> y_i = y.block(i * x_height, 0, x_height, 1);

		Ref<RowMatrixXd> temp_y = Cy.block(height_total * i, 0, height_total, 1);
		temp_y.block(0, 0, Dcut_height, 1) = Dy_i.block(Dcut[0], 0, Dcut_height, 1);
		temp_y.block(Dcut_height, 0, height_total - Dcut_height, 1) =
			y_i.block(cut[0], 0, height_total - Dcut_height, 1);

		if (transformation == "fod")
			temp_y(0, 0) = NAN;
		Ref<RowMatrixXd> Dx_i = Dx.block(i * Dx_height, 0, Dx_height, Dx.cols());
		Ref<RowMatrixXd> x_i = x.block(i * x_height, 0, x_height, x.cols());
		Ref<RowMatrixXd> temp_x =
			Cx.block(height_total * i, 0, height_total, Cx.cols());
		temp_x.block(0, 0, Dcut_height, width) =
			Dx_i.block(Dcut[0], 0, Dcut_height, width);

		temp_x.block(Dcut_height, 0, height_total - Dcut_height, width) =
			x_i.block(cut[0], 0, height_total - Dcut_height, width);

		temp_x.block(0, width, Dcut_height, 1) = RowMatrixXd::Zero(Dcut_height, 1);

		temp_x.block(Dcut_height, width, height_total - Dcut_height, 1) =
			RowMatrixXd::Ones(height_total - Dcut_height, 1);
	}
	return std::make_tuple(Cy, Cx);
}

std::tuple<RowMatrixXd, RowMatrixXd>
get_final_xy_diffGMM(Ref<RowMatrixXd> Dy, Ref<RowMatrixXd> Dx, vector<int> Dcut, int Dcut_height)
{
	int Dx_height = Dx.rows() / N;

	int width = Dx.cols();

	RowMatrixXd Cx(Dcut_height * N, width);
	RowMatrixXd Cy(Dcut_height * N, 1);
#pragma omp parallel for
	for (int i = 0; i < N; ++i)
	{
		Ref<RowMatrixXd> Dy_i = Dy.block(i * Dx_height, 0, Dx_height, 1);
		Ref<RowMatrixXd> temp_y =
			Cy.block(Dcut_height * i, 0, Dcut_height, 1);
		temp_y.block(0, 0, Dcut_height, 1) =
			Dy_i.block(Dcut[0], 0, Dcut_height, 1);
		Ref<RowMatrixXd> Dx_i =
			Dx.block(i * Dx_height, 0, Dx_height, width);

		Ref<RowMatrixXd> temp_x =
			Cx.block(Dcut_height * i, 0, Dcut_height, width);
		temp_x.block(0, 0, Dcut_height, width) =
			Dx_i.block(Dcut[0], 0, Dcut_height, width);
	}
	return std::make_tuple(Cy, Cx);
}

void prepare_reg_fod(Ref<RowMatrixXd> Diff_x, Ref<RowMatrixXd> Diff_y)
{
	// int xy_height = Diff_x.rows() / N;

	vector<bool> row_if_nan_y = row_has_nan(Diff_y);
	vector<bool> row_if_nan_x = row_has_nan(Diff_x);
	for (std::size_t i = 0, max = row_if_nan_y.size(); i != max; ++i)
		if ((row_if_nan_y[i]) || (row_if_nan_x[i]))
		{
			Diff_x.row(i) = RowMatrixXd::Zero(1, Diff_x.cols());
			Diff_y(i, 0) = 0;
		}
}

std::tuple<int, int, int, double>
prepare_reg(Ref<RowMatrixXd> z_list, Ref<RowMatrixXd> Cx, Ref<RowMatrixXd> Cy,
			z_info z_information, string transformation, bool level)
{
	int xy_height = Cy.rows() / N;
	int z_height = z_list.rows() / N;
	int Cx_width = Cx.cols();
	int z_width = z_list.cols();

	int num_NA = 0;
	// int total = 0;
	int max = 0;
	int min = 0;
	if (transformation == "fod")
		prepare_reg_fod(final_xy_tables["Diff_x"], final_xy_tables["Diff_y"]);
	int start_row;//, end_row;
	//end_row = z_information.z_width;
	start_row = 0;
	if (level)
		start_row = z_information.diff_width;


	for (int i = 0; i < N; ++i)
	{
		Ref<RowMatrixXd> x = Cx.block(i * xy_height, 0, xy_height, Cx_width);
		Ref<RowMatrixXd> y = Cy.block(i * xy_height, 0, xy_height, 1);

		Ref<RowMatrixXd> z = z_list.block(i * z_height, 0, z_height, z_width);

		vector<bool> row_if_nan_x = row_has_nan(x);
		vector<bool> row_if_nan_y = row_has_nan(y);
		int temp = 0;
		// for (int row_id = start_row; row_id < end_row; ++row_id)
		// {
		// 	if ((row_if_nan_x[i]) || (row_if_nan_y[i]))
		// 		temp += 1;
		// }
		// num_NA += temp;
		// if (temp > max)
		// 	max = temp;
		// if (temp < min)
		// 	min = temp;

		//    na_list.append(row_if_nan)
		for (int j = 0; j < z_width; ++j)
			if ((row_if_nan_x[j]) || (row_if_nan_y[j]))
			{
				if (j>=start_row)
					temp+=1;
				x.row(j) = RowMatrixXd::Zero(1, Cx_width);
				y(j, 0) = 0;
				z.col(j) = RowMatrixXd::Zero(z_height, 1);
				// std::cout << i << " " << j << std::endl;
			}
		num_NA += temp;
		if (temp > max)
			max = temp;
		if (temp < min)
			min = temp;
	}
	int width;
	if (level)
		width = z_information.level_width;
	else
		width = z_information.diff_width;
	int nobs = width * N - num_NA;
	int max_obs = width - min;
	int min_obs = width - max;
	double avg_obs = width*1.0 - (1.0 / N)*num_NA;
	//std::cout << width << " " << num_NA << " " << min << " " << avg_obs << std::endl;

	return std::make_tuple(nobs, max_obs, min_obs, avg_obs);
}
