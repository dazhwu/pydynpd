import math

import numpy as np
from scipy import stats

from pydynpd.common_functions import lag
from pydynpd.info import hansen_test_info, AR_test_info


def hansen_overid(W2_inv, N, zs, num_instru, num_indep):
    hansen_test = np.linalg.multi_dot([zs.transpose(), W2_inv, zs]) * (1.0 / N)
    df = num_instru - num_indep
    crit = stats.chi2.ppf(q=0.95, df=df)  # if hansen > crit then H0 is not supported
    p_value = 1 - stats.chi2.cdf(hansen_test, df)

    hansen_test = hansen_test_info(float(hansen_test), df, float(p_value), crit)

    return (hansen_test)


def AR_get_diff_XR(model, beta, ori_residual, r0_height):
    N = model.N

    ori_x = model.final_xy_tables['Cx']
    x0_height = int(ori_x.shape[0] / N)

    if model.options.transformation == 'fod':
        diff_y = model.final_xy_tables['Diff_y']
        diff_x = model.final_xy_tables['Diff_x']

        diff_r = diff_y - diff_x @ beta
    elif model.options.level:
        diff_width = model.z_information.diff_width
        num_col = ori_x.shape[1]
        diff_r = np.empty((diff_width * N, 1), dtype=np.float64)
        diff_x = np.empty((diff_width * N, num_col), dtype=np.float64)
        for i in range(N):
            r_i = ori_residual[(i * r0_height):(i * r0_height + r0_height), :]
            diff_r[(i * diff_width):(i * diff_width + diff_width), 0:1] = r_i[0:diff_width, 0:1]

            x_i = ori_x[(i * r0_height):(i * r0_height + r0_height), :]
            diff_x[(i * diff_width):(i * diff_width + diff_width), :] = x_i[0:diff_width, :]
    else:
        diff_x = ori_x
        diff_r = ori_residual

    return (diff_x, diff_r)


def AR_test(model, zs_list, step, m):
    N = model.N
    z_list = model.z_list
    z_height = int(z_list.shape[0] / N)

    current_step = model.step_results[step - 1]

    ori_residual = current_step.residual
    r0_height = int(ori_residual.shape[0] / N)

    M_XZ_W = current_step._M_XZ_W
    vcov = current_step.vcov

    diff_x, diff_r = AR_get_diff_XR(model, current_step.beta, ori_residual, r0_height)
    r_height = int(diff_r.shape[0] / N)
    x_height = int(diff_x.shape[0] / N)
    AR_list = []
    temp = np.zeros((r_height * N, 1), np.float64)
    lag(diff_r, temp, N, 1, 0)
    for j in range(1, m + 1):
        for i in range(N):
            r_i = diff_r[(r_height * i):(r_height * i + r_height), 0:1]
            r_t_i = r_i.transpose()

            lag_res = np.ndarray((r_height, 1), dtype=np.float64)
            lag(r_i, lag_res, 1, j, 0)
            lag_res[np.isnan(lag_res)] = 0
            lag_res_t = lag_res.transpose()

            x = diff_x[(x_height * i):(x_height * i + x_height), :]
            # z = z_list[(z_height * i):(z_height * i + z_height), :]

            # r_whole_i = ori_residual[(i * r0_height):(i * r0_height + r0_height), 0:1]

            d0_temp = lag_res_t @ r_i
            d1_temp = d0_temp @ r_t_i @ lag_res
            EX_temp = lag_res_t @ x

            zs = zs_list[(z_height * i):(z_height * i + z_height), :]

            # temp3_temp = z @ r_whole_i @ r_t_i @ lag_res
            # temp3_temp = zs @ r_t_i @ lag_res
            temp3_temp = zs @ d0_temp.transpose()
            if i == 0:
                d0 = d0_temp
                d1 = d1_temp
                EX = EX_temp
                temp3 = temp3_temp
            else:
                d0 += d0_temp
                d1 += d1_temp
                EX += EX_temp
                temp3 += temp3_temp

        d2 = (-2) * np.linalg.multi_dot([EX, M_XZ_W, temp3])

        d3 = np.linalg.multi_dot([EX, vcov, EX.transpose()])
        try:
            AR_temp = float(d0 / math.sqrt(d1 + d2 + d3))
        except Exception as e:
            raise Exception('AR test failed')

        P_value = stats.norm.sf(abs(AR_temp)) * 2
        new_AR = AR_test_info(j, AR_temp, P_value)
        AR_list.append(new_AR)

    return (AR_list)
