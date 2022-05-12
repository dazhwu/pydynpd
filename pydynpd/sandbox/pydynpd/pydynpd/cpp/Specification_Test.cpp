#include "Specification_Test.h"

Hansen_test_info::Hansen_test_info(double _test, int _df, double _cri, double _p) {
    test_value = _test;
    df = _df;

    critical_value = _cri;
    P_value = _p;
}


AR_test_info::AR_test_info(int _lag, double _AR, double _P) {
    lag = _lag;
    AR = _AR;
    P_value = _P;
}


Hansen_test_info
hansen_overid(const Ref<const RowMatrixXd> &W2_inv, const Ref<const RowMatrixXd> &zs, int num_instru, int num_indep,
              int N) {


    double hansen_test = (zs.transpose() * W2_inv * zs).value() * (1.0 / N);
    int df = num_instru - num_indep;
    boost::math::chi_squared k2(df);
    //double crit=
    //struct Hansen_test_info tbr={hansen_test, df};
    return Hansen_test_info(hansen_test, df,boost::math::pdf(k2, 0.95) , 1- boost::math::cdf(k2,hansen_test));
}


vector<MatrixXd> AR_get_diff_XR(int N, Eigen::Ref <RowMatrixXd> beta, Eigen::Ref <RowMatrixXd> ori_residual,
                                      Eigen::Ref <RowMatrixXd> ori_x, int diff_width,
                                      Eigen::Ref <RowMatrixXd> Diff_y_table, Eigen::Ref <RowMatrixXd> Diff_x_table,
                                      int r0_height, string transformation, bool level) {

    vector <MatrixXd> tbr;
    if (transformation == "fod") {
        MatrixXd diff_r(diff_width * N, 1);
        diff_r = Diff_y_table - Diff_x_table * beta;
        tbr.push_back(Diff_x_table);
        tbr.push_back(diff_r);
        return tbr;
    } else if (level == true) {
        int num_col = ori_x.cols();
        RowMatrixXd diff_r(diff_width * N, 1);
        RowMatrixXd diff_x(diff_width * N, num_col);
        int starting_row = 0, starting_row_original = 0;
        for (int i = 0; i < N; ++i) {
            diff_r.block(starting_row, 0, diff_width, 1) = ori_residual.block(starting_row_original, 0, diff_width, 1);
            diff_x.block(starting_row, 0, diff_width, num_col) = ori_x.block(starting_row_original, 0, diff_width,
                                                                             num_col);
            starting_row += diff_width;
            starting_row_original += r0_height;
        }
        tbr.push_back(diff_x);
        tbr.push_back(diff_r);
        return tbr;
    } else {
        tbr.push_back(ori_x);
        tbr.push_back(ori_residual);
        return tbr;
    }
}

vector<AR_test_info> AR_test(int N, int m, Eigen::Ref <RowMatrixXd> z_list, Eigen::Ref <RowMatrixXd> zs_list,
                                   Eigen::Ref <RowMatrixXd> ori_residual,
                                   Eigen::Ref <RowMatrixXd> M_XZ_W, Eigen::Ref <RowMatrixXd> vcov,
                                   Eigen::Ref <RowMatrixXd> ori_x, Eigen::Ref <RowMatrixXd> beta,
                                   int diff_width, Eigen::Ref <RowMatrixXd> Diff_y_table,
                                   Eigen::Ref <RowMatrixXd> Diff_x_table, string transformation, bool level) {

    int z_height = z_list.rows() / N;

    int r0_height = ori_residual.rows() / N;
    // Eigen::VectorXd AR_list(m);
    std::vector <AR_test_info> AR_list;
    // Eigen::Ref<MatrixXd> diff_x, diff_r;

    vector <MatrixXd> ret = AR_get_diff_XR(N, beta, ori_residual, ori_x, diff_width,
                                           Diff_y_table, Diff_x_table, r0_height, transformation, level);

    RowMatrixXd diff_x = ret[0];
    RowMatrixXd diff_r = ret[1];

    int r_height = diff_r.rows() / N;
    int x_height = diff_x.rows() / N;
    int x_width = diff_x.cols();

    // temp = np.zeros((r_height * N, 1), np.float64)
    // lag(diff_r, temp, N, 1, 0)
    RowMatrixXd EX, temp3;
    double d0, d1, d2, d3;

    // Map<RowMatrixXd> r_i(NULL, r_height,1);
    // Map<RowMatrixXd> x(NULL, x_height, x_width);
    for (int j = 1; j <= m; ++j) {
        for (int i = 0; i < N; ++i) {

            Eigen::Ref <RowMatrixXd> r_i = diff_r.block(r_height * i, 0, r_height, 1);
            RowMatrixXd r_t_i = r_i.transpose();
            RowMatrixXd lag_res(r_height, 1);

            for (int k = 0; k < j; ++k)
                lag_res(k, 0) = 0;
            lag_res.block(j, 0, r_height - j, 1) = r_i.block(0, 0, r_height - j, 1);

            RowMatrixXd lag_res_t = lag_res.transpose();

            Eigen::Ref <RowMatrixXd> x = diff_x.block(x_height * i, 0, x_height, x_width);

            RowMatrixXd d0_temp = lag_res_t * r_i;

            auto d1_temp = d0_temp * r_t_i * lag_res;
            auto EX_temp = lag_res_t * x;

            Ref <RowMatrixXd> zs = zs_list.block(z_height * i, 0, z_height, 1);

            auto temp3_temp = zs * d0_temp.transpose();

            if (i == 0) {
                d0 = d0_temp.value();
                d1 = d1_temp.value();
                EX = EX_temp;
                temp3 = temp3_temp;
            } else {
                d0  += d0_temp.value();
                d1  += d1_temp.value();
                EX.noalias()  += EX_temp;
                temp3.noalias()  += temp3_temp;
            }
        }
        d2 = (-2) * (EX * M_XZ_W * temp3).value();
        d3 = (EX * vcov * EX.transpose()).value();
        double AR_temp = d0 / sqrt(d1 + d2 + d3);
        double P_value = 2 * (1 - standard_normalCDF(abs(AR_temp)));

        AR_test_info new_AR = AR_test_info(j, AR_temp, P_value);
        AR_list.push_back(new_AR);


    }

    return AR_list;
}
