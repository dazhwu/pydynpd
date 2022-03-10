import tokenize


import numpy as np
cimport numpy as np
cimport cython
import pandas as pd
from pandas import DataFrame

import variable

from command import parse_command
from panel_data import new_panel_data
from libc.math cimport pow, sqrt

import time

cdef:
    dict variables
    int N, T
    int diff_width, diff_height, level_width, level_height
    bint level
    np.double_t [:,:]  H, W, W2, XZ, Zy, vcov_step1, vcov_step2
    list z_list, z_t_list
    list Cx_list, Cy_list
    list residual, residual2
    str  method



def pydpd(command_str, df: DataFrame, list identifier):
    cdef:
        dict variables, final_xy_tables
        list z_list

    variables = parse_command(command_str)
    level=True

    z_list, z_information, df_inf, final_xy_tables = new_panel_data(df, ['id', 'year'], variables, level)

    N=df_inf.N

    Cx_list=final_xy_tables['Cx']
    Cy_list=final_xy_tables['Cy']

    prepare_reg(Cx_list, Cy_list, z_list)
    H=get_H(z_list, z_information)


    z_t_list=[z.transpose() for z in z_list]

    W=(1.0/N)*sum_product([z_list, H,  z_t_list], N)


    XZ=sum_product([z_list, Cx_list],N).transpose()

    Zy=sum_product([z_list,Cy_list], N).transpose()


    XZ_W = np.matmul(XZ, np.linalg.pinv(W))

    M = np.linalg.pinv(np.matmul(XZ_W, XZ.transpose()))

    beta = np.linalg.multi_dot([M, XZ_W, Zy.transpose()])
    print('step 1 estimator')
    print(beta)

    residual = calculate_residual(Cy_list, Cx_list, beta)

    residual_t=[mat.transpose() for mat in residual]

    ZuuZ= sum_product([z_list, residual, residual_t, z_t_list], N)
    W2 = ZuuZ * (1.0 / N)


    vcov_step1 = N*np.linalg.multi_dot([M, XZ_W, W2, XZ_W.transpose(),M])

    print('std err')
    for i in range(0, 4):
        print(sqrt(vcov_step1[i, i]))

    XZ_W2 = np.matmul(XZ, np.linalg.pinv(W2))

    M2 = np.linalg.pinv(np.matmul(XZ_W2, XZ.transpose()))

    beta2 = np.linalg.multi_dot([M2, XZ_W2, Zy.transpose()])

    print('step 2 estimator')
    print(beta2)

    residual2 = calculate_residual(Cy_list, Cx_list, beta2)

    zs2=sum_product([z_list, residual2], N)

    vcov_step2=Windmeijer(M2, XZ_W2, W2, zs2, vcov_step1, Cx_list, z_list, residual)

    for i in range(0, 4):
        print(sqrt(vcov_step2[i, i]))

    temp=N*np.linalg.multi_dot([zs2.transpose(), W2, zs2])

    temp=np.linalg.multi_dot([zs2.transpose(),np.linalg.pinv(ZuuZ), zs2])

def calculate_residual(y_list, x_list, beta):
    cdef int i, n_obs=len(y_list)
    
    cdef list tbr=[]

    for i in range(n_obs):
        temp = y_list[i] - np.matmul(x_list[i], beta)
        tbr.append(temp)

    return(tbr)


def sum_product(listOflist, int n_rows):
    cdef int i, j, num_elements=len(listOflist)
    cdef list list_temp, var_list
    
    cdef np.ndarray[np.double_t, ndim = 2] var_mat
    cdef np.ndarray[np.double_t, ndim = 2] temp, tbr

    for i in range(n_rows):
        list_temp=[]
        for j in range(num_elements):
            if type(listOflist[j])==list:
                var_list=listOflist[j]
                list_temp.append(var_list[i])
            elif type(listOflist[j])==np.ndarray:
                var_mat=listOflist[j]
                list_temp.append(var_mat)
            else:
                pass  #throw error
        temp=np.linalg.multi_dot(list_temp)
        if i==0:
            tbr=temp
        else:
            tbr=np.add(tbr, temp)

    return(tbr)

def Windmeijer(np.ndarray[np.double_t, ndim = 2] M2, 
               np.ndarray[np.double_t, ndim = 2] XZ_W2, 
               np.ndarray[np.double_t, ndim = 2] W2, 
               np.ndarray[np.double_t, ndim = 2] zs2,
               np.ndarray[np.double_t, ndim = 2] vcov_step1,
               list Cx_list,
               list z_list,
               list residual):
    
    cdef int i, j, N=len(Cx_list)

    cdef np.ndarray[np.double_t, ndim = 2] D, temp, x,  z, u, xu


    D = np.empty((M2.shape[0], M2.shape[1]), dtype='float64')

    for j in range(0, Cx_list[0].shape[1]):

        for i in range(0, N):
            x = Cx_list[i]

            u = residual[i]
            z = z_list[i]

            xu = np.matmul(x[:, j:(j + 1)], u.transpose())

            temp=np.linalg.multi_dot([z, xu+xu.transpose(), z.transpose()])
            #temp = np.matmul(z, xu + xu.transpose())

            #temp = np.matmul(temp, z.transpose())
            if i == 0:
                zxz = temp
            else:
                zxz += temp

        partial_dir = (-1.0 / N) * zxz

        Dj = np.linalg.multi_dot([M2, XZ_W2, partial_dir, np.linalg.pinv(W2), zs2])
        Dj = (-1) * Dj

        D[:, j:(j + 1)] = Dj


    temp = np.multiply(N, M2) + np.multiply(N,  np.matmul(D, M2))  + np.multiply(N,  np.matmul(M2, D.transpose()))

    temp = temp + np.matmul(np.matmul(D, vcov_step1), D.transpose())

    return(temp)

def get_H(list z_list, z_inf):
    cdef int width
    #cdef numpy.ndarray H

    width = z_list[0].shape[1]

    tbr=np.zeros((width,width), dtype='float64')
    i, j = np.indices(tbr.shape)
    tbr[np.logical_and(i==j,i<z_inf.diff_width)]=2
    tbr[np.logical_and(i==j-1,  j<z_inf.diff_width)]=-1
    tbr[np.logical_and(j==i-1,  i<z_inf.diff_width)]=-1

    tbr[np.logical_and(i == j,  i >= z_inf.diff_width)] = 1

    tbr[np.logical_and(i==j+z_inf.diff_width, j<z_inf.diff_width)]=-1
    tbr[np.logical_and(i == 1+j + z_inf.diff_width, j < z_inf.diff_width)] = 1
    tbr[np.logical_and(j == i + z_inf.diff_width, i < z_inf.diff_width)] = -1
    tbr[np.logical_and(j == 1 + i + z_inf.diff_width, i < z_inf.diff_width)] = 1

    return(tbr)






def prepare_reg(list Cx_list, list Cy_list, list z_list):
    cdef int i, j, N=len(Cy_list)
    cdef list na_list=[]
    cdef np.ndarray[np.double_t, ndim = 2] x, y, z
    #na_list=[]
    for i in range(N):
        x=Cx_list[i]
        y=Cy_list[i]
        z=z_list[i]
        row_if_nan = np.logical_or ( np.isnan(x).any(axis=1) , np.isnan(y).any(axis=1))

        na_list.append(row_if_nan)
        for j in range(0, len(row_if_nan)):
            if row_if_nan[j] == True:
                x[j, :] = 0
                y[j,:]=0
                z[:,j]=0
                #n_obs = n_obs - 1
    #return(Cx_list, Cy_list, z_list)