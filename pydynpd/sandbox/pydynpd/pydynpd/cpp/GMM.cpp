
#include "GMM.h"

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

std::vector <Step_Result> results;
//std::tuple<vector<Step_Result> , Hansen_test_info , vector<AR_test_info>



std::tuple <struct Regression, Hansen_test_info, vector<AR_test_info>, struct basic_info>
regular_process(Eigen::Ref <RowMatrixXd> z_list,
                Eigen::Ref <RowMatrixXd> Cx_list,
                Eigen::Ref <RowMatrixXd> Cy_list,
				struct basic_info model_info,
				struct model_options options
) {
	
	string transformation=options.transformation;
	int steps=options.steps;
	bool level =options.level;
	int N =model_info.N;
	int T=model_info.T;
	int num_instru=model_info.num_instr;
	int num_indep=model_info.num_indep;
	int num_obs=model_info.num_obs;
	int diff_width=model_info.diff_width;
	
	
	RowMatrixXd _XZ, _Zy;
	
	MatrixXd _z_t_list = z_list.transpose();
	results.clear();
	std::tie(_XZ, _Zy) = calculate_basic(z_list, Cx_list, Cy_list, N);
	RowMatrixXd _XZ_t = _XZ.transpose();
	RowMatrixXd _Zy_t = _Zy.transpose();
	GMM(N, T, num_obs, diff_width, z_list, _z_t_list, Cx_list, Cy_list, _XZ, _XZ_t, _Zy, _Zy_t, 1, transformation,
	    level);
	GMM(N, T, num_obs, diff_width, z_list, _z_t_list, Cx_list, Cy_list, _XZ, _XZ_t, _Zy, _Zy_t, 2, transformation,
	    level);
	
	int num_steps=results.size();
	Ref <RowMatrixXd> _zs_list = results[num_steps - 1]._zs_list;

	//std::tie<Hansen_test_info hansen, vector<AR_test_info> ar_list> =perform_test(Cx_list, Cy_list, z_list, _zs_list, num_instru, num_indep, 2, diff_width,  N, transformation,level);
	Hansen_test_info hansen = perform_hansen_test(2, num_instru, num_indep, N);
	RowMatrixXd reg_table=regression_table(results[num_steps-1].beta,results[steps-1].std_err);
	
	vector <AR_test_info> ar_list = perform_AR_test(Cx_list, Cy_list, z_list, _zs_list, 2, diff_width, N,
	                                                transformation, level);

	Regression tbr = Regression(reg_table, results[num_steps-1].vcov, results[num_steps - 1].W,z_list, Cy_list, Cx_list);
	//return std::make_tuple(results[steps-1].beta,results[steps-1].std_err,results[steps-1].vcov,results[steps-1].SS, hansen, ar_list, z_list);
	
	model_info.SS=results[num_steps-1].SS;
	
	return std::make_tuple(tbr, hansen, ar_list, model_info);
	//     self.perform_test(model, 2)
	// else:
	//     self.iterative_GMM(model, _XZ, _XZ_t, _Zy, _Zy_t)
	//     self.perform_test(model, model.options.steps)
}

RowMatrixXd regression_table(RowMatrixXd beta, RowMatrixXd std_err){
	int num_coe=beta.rows();
	RowMatrixXd tbr(num_coe, 4);
	//std::cout << beta << std::endl;
	for (int i=0; i<num_coe; ++i){
		tbr(i,0)=beta(i,0);
		tbr(i,1)=std_err(i,0);
		if (tbr(i,1) != 0){
			tbr(i,2)=tbr(i,0)/tbr(i,1);
			tbr(i,3)=2 * (1 - standard_normalCDF(abs(tbr(i,2))));
			
		}
	}
	//std::cout<<tbr<<std::endl;
 	return tbr;
 }

void GMM(int N, int T, int num_obs, int diff_width, Ref <RowMatrixXd> z_list,
         Ref <MatrixXd> z_list_t,
         Ref <RowMatrixXd> Cx_list,
         Ref <RowMatrixXd> Cy_list, Ref <RowMatrixXd> _XZ,
         Ref <RowMatrixXd> _XZ_t, Ref <RowMatrixXd> _Zy,
         Ref <RowMatrixXd> _Zy_t, int step, string transformation, bool level) {
	RowMatrixXd W;
	if (step == 1) {
		RowMatrixXd H1;
		H1 = get_H1(z_list, diff_width, T, transformation, level);
		//std::cout << H1 << std::endl;
		W = calculate_W(H1, z_list, z_list_t, 1, N);

	} else {
		Step_Result previous_step = results[step - 2];
		W = previous_step.W_next;

	}
	RowMatrixXd W_inv = common_inv(W);
	RowMatrixXd _XZ_W = _XZ * W_inv;
	RowMatrixXd _M_inv = _XZ_W * _XZ_t;
	RowMatrixXd M = common_inv(_M_inv);
	RowMatrixXd _M_XZ_W = M * _XZ_W;
	////std::cout << W << std::endl;
	RowMatrixXd beta = _M_XZ_W * _Zy_t;
	RowMatrixXd residual = calculate_residual(Cy_list, Cx_list, beta);
	RowMatrixXd _residual_t = residual.transpose();
	double SS = (_residual_t * residual).value() * (1.0 / 2 / num_obs);

	int z_height = z_list.rows() / N;
	int z_width = z_list.cols();
	int r_height = residual.rows() / N;
	//int r_width = residual.cols(); //=1
	RowMatrixXd _zs_list(N * z_height, 1);
	RowMatrixXd zs(z_height, 1), ZuuZ(z_height, z_height);
	for (int i = 0; i < N; ++i) {
		//Map<RowMatrixXd> z (z_list.block(i * z_height, 0, z_height, z_width));
		Ref <RowMatrixXd> u = residual.block(i * r_height, 0, r_height, 1);
		_zs_list.block(i * z_height, 0, z_height, 1) = z_list.block(i * z_height, 0, z_height, z_width) * u;
		Ref <RowMatrixXd> temp_zs = _zs_list.block(i * z_height, 0, z_height, 1);
		if (i == 0) {
			zs = temp_zs;
			ZuuZ = temp_zs * temp_zs.transpose();
		} else {
			zs.noalias()  += temp_zs;
			ZuuZ.noalias()  += temp_zs * temp_zs.transpose();
		}
	}
	RowMatrixXd W_next = ZuuZ * (1.0 / N);

	RowMatrixXd _vcov;
	if (step == 1)
		_vcov = vcov_step_1(_M_XZ_W, W_next, N);
	else
		_vcov = vcov(z_list, Cx_list, M, _M_XZ_W, W_inv, zs, step, N);
	int vcov_width = _vcov.rows();
	RowMatrixXd std_err(vcov_width, 1), temp;
	temp = _vcov.diagonal();
	for (int i = 0; i < vcov_width; ++i)
		std_err(i, 0) = sqrt(temp(i, 0));
	
	Step_Result current_step = Step_Result(W, W_inv, W_next, _XZ_W, M, zs, ZuuZ, _vcov, _M_XZ_W, _zs_list, beta,
	                                       std_err, residual, _residual_t, SS);
	
	results.push_back(current_step);
	
}

std::tuple <RowMatrixXd, RowMatrixXd> calculate_basic(Ref <RowMatrixXd> z_list,
                                                      Ref <RowMatrixXd> Cx_list,
                                                      Ref <RowMatrixXd> Cy_list, int N) {
	int z_height = z_list.rows() / N;
	int x_height = Cx_list.rows() / N;
	int x_width = Cx_list.cols();
	int z_width = z_list.cols();
	RowMatrixXd temp_xz, temp_zy;
	RowMatrixXd zx;
	for (int i = 0; i < N; ++i) {
		Ref <RowMatrixXd> z = z_list.block(z_height * i, 0, z_height, z_width);
		//Ref<RowMatrixXd> z_t = _z_t_list.block(0, z_height * i, z_width, z_height);

		Ref <RowMatrixXd> x = Cx_list.block(x_height * i, 0, x_height, x_width);

		Ref <RowMatrixXd> y = Cy_list.block(x_height * i, 0, x_height, 1);
		zx = z * x;
		if (i == 0) {
			temp_xz = zx.transpose();
			temp_zy = (z * y).transpose();
		} else {
			temp_xz.noalias()  += zx.transpose();
			temp_zy.noalias()  += (z * y).transpose();
		}
	}

	return std::make_tuple(temp_xz, temp_zy);
}

RowMatrixXd calculate_W(Ref <RowMatrixXd> H,
                        Ref <RowMatrixXd> z_list,
                        Ref <MatrixXd> _z_t_list, int step,
                        int N) {
	int z_height = z_list.rows() / N;
	int z_width = z_list.cols();
	RowMatrixXd temp_W;
	for (int i = 0; i < N; ++i) {
		Ref <RowMatrixXd> z = z_list.block(z_height * i, 0, z_height, z_width);
		if (i == 0)
			temp_W = z * H * z.transpose();
		else
			temp_W += z * H * z.transpose();
	}
	return temp_W;
}

RowMatrixXd calculate_residual(Ref <RowMatrixXd> y_list,
                               Ref <RowMatrixXd> x_list,
                               Ref <RowMatrixXd> beta) {
	RowMatrixXd tbr = y_list - x_list * beta;
	////std::cout << tbr << std::endl;
	return tbr;
}

RowMatrixXd vcov_step_1(Ref <RowMatrixXd> _M_XZ_W, Ref <RowMatrixXd> W2, int N) {
	//Ref<RowMatrixXd> W2 = step_1.W_next;
	RowMatrixXd tbr = N * (_M_XZ_W * W2 * _M_XZ_W.transpose());
	return tbr;

}

RowMatrixXd vcov(Ref <RowMatrixXd> z_list,
                 Ref <RowMatrixXd> Cx, Ref <RowMatrixXd> M2, Ref <RowMatrixXd> _M2_XZ_W2, Ref <RowMatrixXd> _W2_inv,
                 Ref <RowMatrixXd> zs2,
                 int step, int N) {
	RowMatrixXd tbr;

	Step_Result previous_step = results[step - 2];
	Ref <RowMatrixXd> vcov_step_previous = previous_step.vcov;
	Ref <RowMatrixXd> residual1_t = previous_step._residual_t;

	tbr = Windmeijer(M2, _M2_XZ_W2, _W2_inv, zs2, vcov_step_previous, Cx,
	                 z_list, residual1_t, N);

	return tbr;
}

Hansen_test_info perform_hansen_test(int step, int num_instru, int num_indep, int N) {
	Step_Result step1 = results[0];
	Step_Result step2 = results[1];
	Step_Result current_step = step2;   //to be changed
//	if (step <= 2)
//		current_step = step2;
//	else
//		current_step = results[step - 1];

	Ref <RowMatrixXd> W2_inv = current_step.W_inv;
	Ref <RowMatrixXd> zs = current_step.zs;

	Hansen_test_info hansen = hansen_overid(W2_inv, zs, num_instru, num_indep, N);
	return hansen;

}

vector <AR_test_info> perform_AR_test(Ref <RowMatrixXd> Cx, Ref <RowMatrixXd> Cy,
                                      Ref <RowMatrixXd> z_list, Ref <RowMatrixXd> _zs_list,
                                      int step, int diff_width, int N, string transformation, bool level) {
	Step_Result step2 = results[1];
	Step_Result current_step = step2;   //to be changed
	if (step <= 2)
		current_step = step2;
	else
		current_step = results[step - 1];
	int m = 2;
	Eigen::Ref <RowMatrixXd> Diff_y = Cy;
	Eigen::Ref <RowMatrixXd> Diff_x = Cx;

	vector <AR_test_info> ar_list = AR_test(
			N, m, z_list, _zs_list, current_step.residual,
			current_step._M_XZ_W, current_step.vcov, Cx,
			current_step.beta, diff_width, Diff_y, Diff_x,
			transformation, level);

	return ar_list;
}

RowMatrixXd get_H1(Ref <RowMatrixXd> z_list, int diff_width, int T,
                   string transformation, bool level) {
	int width = z_list.cols();
	RowMatrixXd tbr=RowMatrixXd::Zero(width, width);
	if (transformation == "fd") {
		for (int i = 0; i < diff_width; ++i) {
			tbr(i, i) = 2;
			if (i >= 1)
				tbr(i - 1, i) = -1;
			if (i < diff_width - 1)
				tbr(i + 1, i) = -1;
		}
		if (width > diff_width) {
			for (int i = diff_width; i < width; ++i)
				tbr(i, i) = 1;
			Eigen::Ref <RowMatrixXd> low_left =
					tbr.block(diff_width, 0, width - diff_width, diff_width);
			for (int i = 0; i < diff_width; ++i) {
				low_left(i, i) = -1;
				low_left(i + 1, i) = 1;
			}
			Eigen::Ref <RowMatrixXd> up_right =
					tbr.block(0, diff_width, diff_width, width - diff_width);
			for (int i = 0; i < diff_width; ++i) {
				up_right(i, i) = -1;
				up_right(i, i + 1) = 1;
			}
		}
	} else {
		tbr = get_H1_fod(width, diff_width, T, level);
	}
	return tbr;
}

RowMatrixXd get_H1_fod(int width, int diff_width, int T, bool level) {
	RowMatrixXd D_up = generate_D_matrix(diff_width, T);
	if (level) {
		RowMatrixXd D(width, T);
		D.block(0, 0, diff_width, T) = D_up;
		Eigen::Ref <RowMatrixXd> low_right =
				D.block(diff_width, T - (width - diff_width), width - diff_width, width - diff_width);
		for (int i = 0; i < width - diff_width; ++i) {
			low_right(i, i) = 1;
		}
		return D * D.transpose();
	} else {
		return D_up * D_up.transpose();
	}
}

RowMatrixXd generate_D_matrix(int height, int T) {
	RowMatrixXd D(height, T);
	for (int i = 0; i < height; ++i) {
		for (int j = i; j < T; ++j) {
			if (i == j)
				D(i, j) = sqrt((T - i - 1) / (T - i));
			else
				D(i, j) = (-1) * sqrt(1 / ((T - i) * (T - i - 1)));
		}
	}
	return D;
}

RowMatrixXd Windmeijer(Ref <RowMatrixXd> M2, Ref <RowMatrixXd> _M2_XZ_W2,
                       Ref <RowMatrixXd> W2_inv, Ref <RowMatrixXd> zs2,
                       Ref <RowMatrixXd> vcov_step1, Ref <RowMatrixXd> Cx_list,
                       Ref <RowMatrixXd> z_list, Ref <RowMatrixXd> residual1_t, int N) {
	RowMatrixXd D(M2.rows(), M2.cols());

	int x_height = Cx_list.rows() / N;
	int z_height = z_list.rows() / N;
	int x_width = Cx_list.cols();
	int z_width = z_list.cols();

	RowMatrixXd zxz;
	//std::cout << "4.1" << std::endl;
	for (int j = 0; j < x_width; ++j) {
		for (int i = 0; i < N; ++i) {
			Ref <RowMatrixXd> x = Cx_list.block(i * x_height, j, x_height, 1);
			//std::cout << residual1_t.cols() << std::endl;
			Ref <RowMatrixXd> u_t = residual1_t.block(0, i * x_height, 1, x_height);

			//std::cout << "4.2" << std::endl;
			Ref <RowMatrixXd> z = z_list.block(i * z_height, 0, z_height, z_width);

			RowMatrixXd xu = x * u_t;
			//RowMatrixXd temp = z * (xu + xu.transpose()) * z.transpose();
			if (i == 0)
				zxz = z * (xu + xu.transpose()) * z.transpose();
			else
				zxz.noalias()  += z * (xu + xu.transpose()) * z.transpose();
		}
		RowMatrixXd partial_dir = (-1.0 / N) * zxz;

		RowMatrixXd Dj = _M2_XZ_W2 * partial_dir * W2_inv * zs2;
		Dj = (-1) * Dj;
		D.col(j) = Dj;

	}
	RowMatrixXd temp_D_M2 = D * M2;
	RowMatrixXd temp = N * M2 + N * temp_D_M2 + N * temp_D_M2.transpose();
	temp += D * vcov_step1 * D.transpose();
	return temp;

}

