import math

import numpy as np
from scipy.sparse import csc_matrix


def lag(mat, lagged, N, lag_number, fill=np.NaN):
    height = int(mat.shape[0] / N)
    for i in range(N):
        start_row = i * height
        end_row = start_row + height
        mat_i = mat[start_row:end_row, :]
        lagged_i = lagged[start_row:end_row, :]

        lagged_i[0:lag_number, :] = fill
        lagged_i[lag_number:height, :] = mat_i[0:(height - lag_number), :]


def get_first_diff_table(ori_arr: np.ndarray, N: int):
    num_cols = ori_arr.shape[1]
    num_rows = ori_arr.shape[0]
    height = int(num_rows / N)

    lag_arr = np.zeros((num_rows, num_cols), dtype='float64')
    tbr_arr = np.zeros((num_rows, num_cols), dtype='float64')

    lag(ori_arr, lag_arr, N, 1)

    tbr_arr = ori_arr - lag_arr
    return tbr_arr


def get_fod_table(ori_arr: np.ndarray, N: int):
    num_rows = ori_arr.shape[0]
    height = int(num_rows / N)

    num_cols = ori_arr.shape[1]

    tbr = np.empty((num_rows, num_cols), dtype='float64')
    next_sum = np.empty((1, num_cols), dtype='float64')
    this_sum = np.empty((1, num_cols), dtype='float64')
    this_avg = np.empty((1, num_cols), dtype='float64')
    temp = np.empty((height, num_cols), dtype='float64')

    tbr[:] = np.NaN

    this_sum[:] = np.NaN

    for i in range(N):
        ori_i = ori_arr[i * height:(i * height + height), :]
        tbr_i = tbr[i * height:(i * height + height), :]
        temp.fill(np.NaN)
        next_sum.fill(np.NaN)
        next_count = 0
        for j in range(height - 2, -1, -1):

            if np.isnan(ori_i[range(j + 1, j + 2), :]).any(axis=1):
                this_count = next_count
                this_sum = next_sum
                temp[j, :] = temp[j + 1, :]
            else:
                this_count = next_count + 1

                this_sum = np.nansum(np.vstack([next_sum, ori_i[j + 1, :]]), axis=0)
                this_avg = this_sum * (1.0 / this_count)
                temp[j, :] = (ori_i[j, :] - this_avg) * math.sqrt(this_count / (this_count + 1))

            next_sum = this_sum
            next_count = this_count

        tbr_i[0, :] = np.NaN
        tbr_i[range(1, height), :] = temp[range(0, height - 1), :]

    return tbr


def sum_product(listOflist, n_rows):
    num_elements = len(listOflist)

    for i in range(n_rows):
        list_temp = []
        for j in range(num_elements):
            if type(listOflist[j]) == list:
                var_list = listOflist[j]
                list_temp.append(var_list[i])
            elif type(listOflist[j]) == np.ndarray:
                var_mat = listOflist[j]
                list_temp.append(var_mat)
            else:
                pass  # throw error
        temp = np.linalg.multi_dot(list_temp)
        if i == 0:
            tbr = temp
        else:
            tbr += temp

    return (tbr)


def Windmeijer(M2, _M2_XZ_W2, W2_inv, zs2, vcov_step1, Cx_list, z_list, residual1, N):
    D = np.empty((M2.shape[0], M2.shape[1]), dtype='float64')

    x_height = int(Cx_list.shape[0] / N)
    z_height = int(z_list.shape[0] / N)
    for j in range(0, Cx_list.shape[1]):

        for i in range(0, N):
            x = Cx_list[(i * x_height):(i * x_height + x_height), :]

            u = residual1[(i * x_height):(i * x_height + x_height), 0:1]
            z = z_list[(i * z_height):(i * z_height + z_height), :]

            xu = np.matmul(x[:, j:(j + 1)], u.transpose())
            temp = z @ (xu + xu.transpose()) @ z.transpose()
            # temp_zxuzt=z@ xu @ z.transpose()
            # temp=temp_zxuzt + temp_zxuzt.transpose()

            if i == 0:
                zxz = temp
            else:
                zxz += temp

        partial_dir = (-1.0 / N) * zxz

        Dj = np.linalg.multi_dot([_M2_XZ_W2, partial_dir, W2_inv, zs2])
        Dj = (-1) * Dj

        D[:, j:(j + 1)] = Dj

    # temp = np.multiply(N, M2) + np.multiply(N, np.matmul(D, M2)) + np.multiply(N, np.matmul(M2, D.transpose()))
    temp_D_M2 = D @ M2
    temp = np.multiply(N, M2) + np.multiply(N, temp_D_M2) + np.multiply(N, temp_D_M2.transpose())
    temp = temp + np.matmul(np.matmul(D, vcov_step1), D.transpose())
    #
    return (temp)


def make_sparse_list(arr_list):
    nrow = len(arr_list)
    new_list = []
    for i in range(nrow):
        arr = arr_list[i]
        new_arr = csc_matrix(arr)
        new_list.append(new_arr)

    return (new_list)
