import tokenize



import numpy as np
cimport numpy as np
cimport cython
import pandas as pd
from pandas import DataFrame

from variable import gmm_var,  variable
from command import command
from panel_data import panel_data
import time




class dpd:
    def __init__(self, command_str, df: DataFrame, identifier):  # constructor method

        my_command = command(command_str)
        self.variables=my_command.variables
        self.level=True

        data = panel_data(df, ['id', 'year'], my_command.variables, self.level)
        self.__dict__.update(data.__dict__)
        self.data=data

        self.Num_xy=len(self.variables["dep_indep"])
        self.prepare_reg()

        self.H=self.get_H()

        self.z_t_list=[z.transpose() for z in self.z_list]

        self.W=(1.0/self.N)*self.sum_product([self.z_list, self.H,  self.z_t_list], self.N)

        self.XZ=self.sum_product([self.z_list, self.Cx_list],self.N).transpose()

        self.Zy=self.sum_product([self.z_list,self.Cy_list], self.N).transpose()


        XZ_W = np.matmul(self.XZ, np.linalg.pinv(self.W))
        self.M = np.linalg.pinv(np.matmul(XZ_W, self.XZ.transpose()))

        self.beta = np.linalg.multi_dot([self.M, XZ_W, self.Zy.transpose()])
        print('step 1 estimator')
        print(self.beta)

        self.residual = self.calculate_residual(self.Cy_list, self.Cx_list, self.beta)

        self.residual_t=[mat.transpose() for mat in self.residual]

        ZuuZ=self.sum_product([self.z_list, self.residual, self.residual_t, self.z_t_list], self.N)
        self.W2 = ZuuZ * (1.0 / self.N)


        self.vcov_step1 = self.N*np.linalg.multi_dot([self.M, XZ_W, self.W2, XZ_W.transpose(),self.M])

        print('std err')
        for i in range(0, 4):
            print(np.sqrt(self.vcov_step1[i, i]))

        self.XZ_W2 = np.matmul(self.XZ, np.linalg.pinv(self.W2))

        self.M2 = np.linalg.pinv(np.matmul(self.XZ_W2, self.XZ.transpose()))

        self.beta2 = np.linalg.multi_dot([self.M2, self.XZ_W2, self.Zy.transpose()])

        print('step 2 estimator')
        print(self.beta2)

        self.residual2 = self.calculate_residual(self.Cy_list, self.Cx_list, self.beta2)

        self.zs2=self.sum_product([self.z_list, self.residual2], self.N)

        self.vcov_step2=self.Windmeijer(self.M2, self.XZ_W2, self.W2, self.zs2)

        for i in range(0, 4):
            print(np.sqrt(self.vcov_step2[i, i]))
            
        temp=self.N*np.linalg.multi_dot([self.zs2.transpose(), self.W2, self.zs2])

        temp=np.linalg.multi_dot([self.zs2.transpose(),np.linalg.pinv(ZuuZ), self.zs2])

    def calculate_residual(self, y_list, x_list, beta):
        cdef int N
        N=len(y_list)
        tbr=[]

        for i in range(0, N):
            temp = y_list[i] - np.matmul(x_list[i], beta)
            tbr.append(temp)

        return(tbr)


    def sum_product(self, listOflist, int N):
        cdef int i, j, num_elements=len(listOflist)
        cdef list list_temp, var_list
        cdef np.ndarray[np.double_t, ndim = 2] var_mat
        cdef np.ndarray[np.double_t, ndim = 2] temp, tbr

        for i in range(N):
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

    def Windmeijer(self, np.ndarray[np.double_t, ndim = 2] M2, np.ndarray[np.double_t, ndim = 2] XZ_W2, np.ndarray[np.double_t, ndim = 2] W2, np.ndarray[np.double_t, ndim = 2] zs2):
        cdef int i, j

        cdef np.ndarray[np.double_t, ndim = 2] D, temp, x,  z, u, xu


        D = np.empty((M2.shape[0], M2.shape[1]), dtype='float64')

        for j in range(0, self.Cx_list[0].shape[1]):

            for i in range(0, self.N):
                x = self.Cx_list[i]

                u = self.residual[i]
                z = self.z_list[i]

                xu = np.matmul(x[:, j:(j + 1)], u.transpose())

                temp=np.linalg.multi_dot([z, xu+xu.transpose(), z.transpose()])
                #temp = np.matmul(z, xu + xu.transpose())

                #temp = np.matmul(temp, z.transpose())
                if i == 0:
                    zxz = temp
                else:
                    zxz += temp

            partial_dir = (-1.0 / self.N) * zxz

            Dj = np.linalg.multi_dot([M2, XZ_W2, partial_dir, np.linalg.pinv(W2), zs2])
            Dj = (-1) * Dj

            D[:, j:(j + 1)] = Dj


        temp = np.multiply(self.N, M2) + np.multiply(self.N,  np.matmul(D, M2))  + np.multiply(self.N,  np.matmul(M2, D.transpose()))

        temp = temp + np.matmul(np.matmul(D, self.vcov_step1), D.transpose())

        return(temp)

    def get_H(self):
        cdef int width
        #cdef numpy.ndarray H

        width = self.z_list[0].shape[1]

        H=np.zeros((width,width), dtype='float64')
        i, j = np.indices(H.shape)
        H[np.logical_and(i==j,i<self.data.diff_width)]=2
        H[np.logical_and(i==j-1,  j<self.data.diff_width)]=-1
        H[np.logical_and(j==i-1,  i<self.data.diff_width)]=-1

        H[np.logical_and(i == j,  i >= self.data.diff_width)] = 1

        H[np.logical_and(i==j+self.data.diff_width, j<self.data.diff_width)]=-1
        H[np.logical_and(i == 1+j + self.data.diff_width, j < self.data.diff_width)] = 1
        H[np.logical_and(j == i + self.data.diff_width, i < self.data.diff_width)] = -1
        H[np.logical_and(j == 1 + i + self.data.diff_width, i < self.data.diff_width)] = 1

        return(H)


        # V=np.zeros((self.data.diff_width+self.data.level_width, self.data.level_width), dtype='float64')
        # i, j = np.indices(V.shape)
        # V[np.logical_and(i==j, i<self.data.diff_width)]=-1
        # V[np.logical_and(i == j-1, i < self.data.diff_width)] = 1
        # V[i==j+self.data.diff_width]=1
        # V_list=[]
        # for i in range(self.N):
        #     V_list.append(V)
        # return (V_list)
        #H2=np.linalg.multi_dot([V, V.transpose()])
        # if ((H==H2).all()):
        #     print('-------+++++++++++++-------------------------------')



    def prepare_reg(self):
        cdef i, j
        cdef list na_list=[]
        cdef np.ndarray[np.double_t, ndim = 2] x, y, z
        #na_list=[]
        for i in range(0, self.N):
            x=self.Cx_list[i]
            y=self.Cy_list[i]
            z=self.z_list[i]
            row_if_nan = np.logical_or ( np.isnan(x).any(axis=1) , np.isnan(y).any(axis=1))

            na_list.append(row_if_nan)
            for j in range(0, len(row_if_nan)):
                if row_if_nan[j] == True:
                    x[j, :] = 0
                    y[j,:]=0
                    z[:,j]=0
                    #self.n_obs = self.n_obs - 1