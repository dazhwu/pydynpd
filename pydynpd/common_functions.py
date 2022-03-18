
import numpy as np
from scipy.sparse import csc_matrix, csr_matrix
import sys


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

def Windmeijer(M2, _M2_XZ_W2, W2_inv, zs2, vcov_step1, Cx_list, z_list, residual1):
    N = len(Cx_list)

    D = np.empty((M2.shape[0], M2.shape[1]), dtype='float64')

    for j in range(0, Cx_list[0].shape[1]):

        for i in range(0, N):
            x = Cx_list[i]

            u = residual1[i]
            z = z_list[i]

            xu = np.matmul(x[:, j:(j + 1)], u.transpose())

            temp = np.linalg.multi_dot([z, xu + xu.transpose(), z.transpose()])
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

def sum_product_sparse(listOflist, n_rows):
# first element on listOflist must be a list


    for i in range(n_rows):
        mat=listOflist[0][i]

        for j in range(1, len(listOflist)):

            element=listOflist[j]

            if type(element)==list:
                current_mat=element[i]

            elif type(element)==np.ndarray:

                current_mat = element
            else:
                print('not accepted type')
                sys.exit()


            temp=csc_matrix.dot(mat, current_mat)

            mat=temp

        if i==0:
            tbr=temp
        else:
            tbr+=temp
    if not type(tbr)==np.ndarray:
        tbr=tbr.todense()

    return tbr


def make_sparse_list(arr_list):
    nrow=len(arr_list)
    new_list=[]
    for i in range(nrow):
        arr=arr_list[i]
        new_arr=csc_matrix(arr)
        new_list.append(new_arr)

    return(new_list)


