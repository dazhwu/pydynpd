
import numpy as np
from scipy.sparse import csc_matrix, csr_matrix
import sys
from pydynpd.info import sumproduct_task

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

def Windmeijer(M2, _M2_XZ_W2, W2_inv, zs2, vcov_step1, Cx_list, z_list, residual1,N):


    D = np.empty((M2.shape[0], M2.shape[1]), dtype='float64')

    x_height=int(Cx_list.shape[0]/N)
    z_height=int(z_list.shape[0]/N)
    for j in range(0, Cx_list.shape[1]):

        for i in range(0, N):
            x = Cx_list[(i*x_height):(i*x_height+x_height),:]

            u = residual1[(i*x_height):(i*x_height+x_height),0:1]
            z = z_list[(i*z_height):(i*z_height+z_height),:]

            xu = np.matmul(x[:, j:(j + 1)], u.transpose())

            temp = z @ (xu + xu.transpose()) @ z.transpose()
            # temp = np.matmul(z, xu + xu.transpose())

            # temp = np.matmul(temp, z.transpose())
            if i == 0:
                zxz = temp
            else:
                zxz += temp

        partial_dir = (-1.0 / N) * zxz

        Dj = np.linalg.multi_dot([_M2_XZ_W2, partial_dir, W2_inv, zs2])
        Dj = (-1) * Dj

        D[:, j:(j + 1)] = Dj

    temp = np.multiply(N, M2) + np.multiply(N, np.matmul(D, M2)) + np.multiply(N, np.matmul(M2, D.transpose()))

    temp = temp + np.matmul(np.matmul(D, vcov_step1), D.transpose())
    #
    return (temp)



def make_sparse_list(arr_list):
    nrow=len(arr_list)
    new_list=[]
    for i in range(nrow):
        arr=arr_list[i]
        new_arr=csc_matrix(arr)
        new_list.append(new_arr)

    return(new_list)


