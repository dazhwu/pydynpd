import numpy as np
from numpy import ndarray
from numpy.linalg import pinv, multi_dot
import scipy.stats as st
from pandas import DataFrame
from prettytable import PrettyTable

from pydynpd.command import parse_command
from pydynpd.panel_data import new_panel_data
import pydynpd.specification_tests as tests
import pydynpd.common_functions
from pydynpd.common_functions import sum_product, Windmeijer
import time


class abond:

    def __init__(self, command_str, df: DataFrame, identifier: list):

        self.initiate_properties()

        self.variables, options = parse_command(command_str)


        self.twosteps = options.twosteps
        self.level = options.level
        self.timedumm=options.timedumm
        self.collapse=options.collapse

        self.identifier=identifier
        robust = True

        self.z_list, self.z_information, df_inf, final_xy_tables \
            = new_panel_data(df, ['id', 'year'], self.variables, options)

        self._z_t_list = [z.transpose() for z in self.z_list]
        self.num_instru=self.z_information.num_instr

        self.N = df_inf.N

        self.Cx_list = final_xy_tables['Cx']
        self.Cy_list = final_xy_tables['Cy']

        self.num_obs = self.prepare_reg()
        self.H1 = self.get_H1()

        # if self.twosteps:
        self.step_1()
        self.step_2()
        # else:
        #     self.step_1()

        self.generate_summary()
    def initiate_properties(self):

        self.z_information=None
        self.AR_list=None
        self.Cx_list=None
        self.Cy_list=None
        self.residual1=None
        self.residual2=None


        self.H1=None
        self.H2=None
        self.M1=None
        self.M2=None
        self.SS1=None
        self.SS2=None
        self.W1=None
        self.W2=None
        self.XZ=None
        self.XZ_W1=None
        self.XZ_W2=None
        self.ZuuZ=None
        self.Zy=None
        self.vcov_step1=None
        self.vcov_step2=None
        self.zs1=None
        self.zs2=None



        self.N=0
        self.num_obs=0
        self.num_instru=0


    def step_1(self):
        N = self.N
        num_obs = self.num_obs
        z_list = self.z_list
        _z_t_list = self._z_t_list
        Cx_list = self.Cx_list
        Cy_list = self.Cy_list
        H1 = self.H1


        W1 = (1.0 / N) * sum_product([z_list, H1, _z_t_list], N)
        XZ = sum_product([z_list, Cx_list], N).transpose()
        Zy = sum_product([z_list, Cy_list], N).transpose()

        XZ_W1 = np.matmul(XZ, np.linalg.pinv(W1))

        M1 = pinv(np.matmul(XZ_W1, XZ.transpose()))

        beta1 = multi_dot([M1, XZ_W1, Zy.transpose()])

        residual1 = self.calculate_residual(Cy_list, Cx_list, beta1)

        _residual1_t = [mat.transpose() for mat in residual1]
        SS1 = sum_product([_residual1_t, residual1], N) / 2 / num_obs
        zs1 = sum_product([z_list, residual1], N)

        ZuuZ = sum_product([z_list, residual1, _residual1_t, _z_t_list], N)
        W2 = ZuuZ * (1.0 / N)

        self.W1 = W1
        self.XZ = XZ
        self.Zy = Zy
        self.XZ_W1 = XZ_W1
        self.M1 = M1
        self.beta1 = beta1
        self.residual1 = residual1
        self._residual1_t = _residual1_t
        self.SS1 = SS1
        self.zs1 = zs1
        self.ZuuZ=ZuuZ
        self.W2=W2

        self.vcov_step1=self.vcov(1)
        self.std_err1 = np.sqrt(np.diag(self.vcov_step1))

    def step_2(self):
        z_list = self.z_list
        _z_t_list = self._z_t_list
        Cx_list = self.Cx_list
        Cy_list = self.Cy_list
        residual1 = self.residual1
        _residual1_t = self._residual1_t
        num_obs = self.num_obs
        XZ = self.XZ
        Zy = self.Zy
        N = self.N
        ZuuZ=self.ZuuZ
        W2=self.W2


        XZ_W2: ndarray = np.matmul(XZ, pinv(W2))

        M2 = pinv(np.matmul(XZ_W2, XZ.transpose()))



        beta2 = multi_dot([M2, XZ_W2, Zy.transpose()])

        residual2 = self.calculate_residual(Cy_list, Cx_list, beta2)
        _residual2_t = [mat.transpose() for mat in residual2]

        SS2 = sum_product([_residual2_t, residual2], N) / 2 / num_obs
        zs2 = sum_product([z_list, residual2], N)


        self.H2 = [np.matmul(r, r.transpose()) for r in self.residual1]
        self.XZ_W2 = XZ_W2
        self.M2 = M2

        self.beta2 = beta2
        self.residual2 = residual2
        self._residual2_t = _residual2_t
        self.SS2 = SS2
        self.zs2 = zs2

        self.vcov_step2 = self.vcov(2)
        self.std_err2 = np.sqrt(np.diag(self.vcov_step2))

    def calculate_residual(self, y_list, x_list, beta):

        tbr = []

        for i in range(self.N):
            temp = y_list[i] - np.matmul(x_list[i], beta)
            tbr.append(temp)

        return (tbr)

    def vcov(self, step: int):
        # report robust vcov only
        if step == 2:
            return Windmeijer(self.M2, self.XZ_W2, self.W2, self.zs2,
                                               self.vcov_step1, self.Cx_list, self.z_list, self.residual1)
        else:
            return self.N * np.linalg.multi_dot([self.M1, self.XZ_W1, self.W2,
                                                 self.XZ_W1.transpose(), self.M1])

        # elif self.twosteps:   #two steps non robust
        #     return(self.M2*self.N)
        # else:                    #one step non robust
        #     ss=sum_product([self.residual1_t,self.residual1], self.N) #*(1/(self.num_obs-self.beta1.shape[0]))
        #     print('-------------')
        #     #print(ss)
        #     return(self.M1*self.N)
        #     #return((self.N)*np.linalg.multi_dot([self.XZ_W1, self.W2, self.XZ.transpose(), self.M1, self.M1]))

    def get_H1(self):
        z_list = self.z_list
        z_inf = self.z_information
        width = z_list[0].shape[1]

        tbr = np.zeros((width, width), dtype='float64')
        i, j = np.indices(tbr.shape)
        tbr[np.logical_and(i == j, i < z_inf.diff_width)] = 2
        tbr[np.logical_and(i == j - 1, j < z_inf.diff_width)] = -1
        tbr[np.logical_and(j == i - 1, i < z_inf.diff_width)] = -1

        tbr[np.logical_and(i == j, i >= z_inf.diff_width)] = 1

        tbr[np.logical_and(i == j + z_inf.diff_width, j < z_inf.diff_width)] = -1
        tbr[np.logical_and(i == 1 + j + z_inf.diff_width, j < z_inf.diff_width)] = 1
        tbr[np.logical_and(j == i + z_inf.diff_width, i < z_inf.diff_width)] = -1
        tbr[np.logical_and(j == 1 + i + z_inf.diff_width, i < z_inf.diff_width)] = 1

        return (tbr)

    def prepare_reg(self):

        Cy_list = self.Cy_list
        Cx_list = self.Cx_list
        z_list = self.z_list

        N = self.N
        na_list = []
        num_NA=0
        #num_NA_diff = 0
        #num_NA_level=0


        # max_observations_per_group = int((np.size(Cy_list[0]) + 1) / 2)

        for i in range(N):
            x = Cx_list[i]
            y = Cy_list[i]
            z = z_list[i]
            row_if_nan = np.logical_or(np.isnan(x).any(axis=1), np.isnan(y).any(axis=1))
            # num_NA+=np.count_nonzero(row_if_nan[range(0,int((np.size(y)+1)/2))])
            # num_NA2 += np.count_nonzero(row_if_nan[range(max_observations_per_group,??? )])
            num_NA += np.count_nonzero(row_if_nan)
            if self.level:
                num_NA -= np.count_nonzero(row_if_nan[range(0, self.z_information.diff_width)])

            na_list.append(row_if_nan)
            for j in range(0, len(row_if_nan)):
                if row_if_nan[j] == True:
                    x[j, :] = 0
                    y[j, :] = 0
                    z[:, j] = 0
                    # n_obs = n_obs - 1
        # return(Cx_list, Cy_list, z_list)
        if self.level:
            return (self.z_information.level_width * N - num_NA)
        else:
            return (self.z_information.diff_width*N-num_NA)


    def generate_summary(self):

        self.hansen = tests.hansen_overid(self.ZuuZ, self.zs2, self.num_instru,
                                          self.Cx_list[0].shape[1])
        self.AR_list = tests.AR_test(self, 2)

        if self.twosteps:
            str_steps='two-step '
        else:
            str_steps='one-step '

        if self.level:
            str_gmm='system GMM'
        else:
            str_gmm='difference GMM'

        print('Dynamic panel-data estimation, ' + str_steps + str_gmm)
        print(self.basic_information())
        print(self.regression_table())
        print(self.test_results())


    def basic_information(self):
        basic_table = PrettyTable()
        basic_table.field_names = ["    ", "   ", "  "]
        basic_table.border=False
        basic_table.header=False
        basic_table.align = 'l'
        basic_table.add_row(['Group variable: ' + self.identifier[0], ' ' , 'Number of obs = ' + str(self.num_obs) ])
        basic_table.add_row(['Time variable: ' + self.identifier[1], ' ', 'Number of groups = ' + str(self.N)])
        basic_table.add_row(['Number of instruments = ' + str(self.num_instru), ' ', ''])

        return(basic_table.get_string())

    def test_results(self):


        str_toprint='Hansen test of overid. restrictions: chi(' +str(self.hansen.df) + ') = ' + '{:.3f}'.format(self.hansen.test_value)
        str_toprint=str_toprint + ' Prob > Chi2 = ' + '{:.3f}'.format(self.hansen.p_value) +'\n'


        for i in range(len(self.AR_list)):

            AR=self.AR_list[i]
            P=st.norm.sf(abs(AR)) * 2
            str_toprint=str_toprint + 'Arellano-Bond test for AR(' + str(i+1) + ') in first differences: z = ' + "{:.2f}".format(AR) + ' Pr > z =' + '{:.3f}'.format(P) + '\n'
        return (str_toprint)

    def regression_table(self):
        if self.twosteps:
            beta=self.beta2
            std_err=self.std_err2
        else:
            beta = self.beta1
            std_err = self.std_err1

        dep_name=self.variables['dep_indep'][0].name
        var_names=[]
        for i in range(1, len(self.variables['dep_indep'])):
            var_name = self.variables['dep_indep'][i].name
            var_lag = self.variables['dep_indep'][i].lag
            if (var_lag) >= 1:
                var_name = 'L' + str(var_lag) + '.' + var_name
            var_names.append(var_name)

        if self.level:
            var_names.append('_con')

        num_indep=len(var_names)
        r_table = PrettyTable()

        r_table.field_names=[dep_name, "coef.", "Corrected Std. Err.", "z", "P>|z|"]

        r_table.float_format = '.7'
        #, "z", "P>|z|", "[95% Conf. Interval]" ]
        for i in range(num_indep):
            var_name=var_names[i]
            coeff=beta[i,0]
            stderr=std_err[i]
            z=coeff/stderr
            p=st.norm.sf(abs(z))*2
            r_table.add_row([var_name, coeff, stderr, z, p])


        return r_table.get_string()