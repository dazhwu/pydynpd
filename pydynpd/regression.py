import numpy as np
from numpy import ndarray
from numpy import matmul
from numpy.linalg import pinv, multi_dot
import scipy.stats as st
from pandas import DataFrame
from prettytable import PrettyTable

from pydynpd.command import parse_command
from pydynpd.panel_data import new_panel_data
import pydynpd.specification_tests as tests
import pydynpd.common_functions
from pydynpd.info import regression_info, sumproduct_task
from pydynpd.common_functions import sum_product, Windmeijer
import time


class abond:

    def __init__(self, command_str, df: DataFrame, identifier: list):

        self.initiate_properties()

        self.variables, options = parse_command(command_str, df.columns)


        self.steps = options.steps

        self.level = options.level
        self.timedumm=options.timedumm
        self.collapse=options.collapse

        self.identifier=identifier
        robust = True

        self.z_list, self.z_information, df_inf, final_xy_tables \
            = new_panel_data(df, identifier, self.variables, options)

        self._z_t_list = self.z_list.transpose()
        self.num_instru=self.z_information.num_instr

        self.N = df_inf.N

        self.Cx_list = final_xy_tables['Cx']
        self.Cy_list = final_xy_tables['Cy']

        self.num_obs = self.prepare_reg()
        
        self.result_list=[]

        self.H1 = self.get_H1()
        self._XZ, self._Zy = self.calculate_basic()
        self._XZ_t=self._XZ.transpose()
        self._Zy_t=self._Zy.transpose()

        self.GMM(1)
        if self.steps==1 or self.steps==2:
              #step 1
            self.GMM(2)  # step 1
        else:
            current_step=1
            converge=False
            while not converge:
                previous_step=current_step
                current_step+=1
                print('step' + str(current_step))
                self.GMM(current_step)
                beta_current=self.result_list[current_step-1].beta
                beta_previous=self.result_list[current_step-2].beta
                for j in range(beta_current.shape[0]):
                    temp=(beta_current[j]-beta_previous[j])**2
                    temp2=(beta_previous[j])**2
                    if j==0:
                        nom=temp
                        denom=temp2
                    else:
                        nom+=temp
                        denom+=temp2
                crit=np.sqrt(nom/denom)

                if crit < 0.00001:
                    converge=True
                    self.steps=current_step
                    print('converged')


        self.generate_summary(self.steps)


    def initiate_properties(self):

        self.N=0
        self.num_obs=0
        self.num_instru=0

    def GMM(self, step):
        N = self.N
        num_obs = self.num_obs
        z_list = self.z_list
        _z_t_list = self._z_t_list
        Cx_list = self.Cx_list
        Cy_list = self.Cy_list
        _XZ=self._XZ
        _XZ_t=self._XZ_t
        _Zy=self._Zy
        _Zy_t=self._Zy_t
        
        if step==1:            
            W = self.calculate_W(self.H1)
            current_step = regression_info(W)
            self.result_list.append(current_step)

        if step>=2:
            current_step = self.result_list[step - 1]
            W = current_step.W

        W_inv = current_step.W_inv
        _XZ_W = _XZ @ W_inv
        
        _M_inv=_XZ_W @ _XZ_t
        M = pinv(_M_inv)
        
        _M_XZ_W = M @ _XZ_W
        
        beta=_M_XZ_W @ _Zy_t

        
        residual = self.calculate_residual(Cy_list, Cx_list, beta)

        _residual_t = residual.transpose()

        SS = (_residual_t @ residual) * (1.0 / 2 / num_obs)

        z_height=int(z_list.shape[0]/N)
        r_height=int(residual.shape[0]/N)
        for i in range(N):
            z=z_list[(i*z_height):(i*z_height+z_height),:]
            u=residual[(i*r_height):(i*r_height+r_height),:]
            #u_t=_residual_t[:, (i*r_height):(i*r_height+r_height)]
            if i==0:
                zs= z @ u
                ZuuZ= zs @ zs.transpose()
            else:
                temp_zs= z @ u
                zs += temp_zs
                ZuuZ += temp_zs @ temp_zs.transpose()


        W_next = ZuuZ * (1.0 / N)

        current_step._XZ_W = _XZ_W
        current_step.M = M
        current_step._M_XZ_W = _M_XZ_W
        current_step.beta = beta
        current_step.residual = residual
        current_step._residual_t = _residual_t
        current_step.SS = SS
        current_step.zs = zs

        current_step.ZuuZ = ZuuZ

        next_step = regression_info(W_next)
        self.result_list.append(next_step)

        current_step.vcov =self.vcov(step)
        current_step.std_err = np.sqrt(np.diag(current_step.vcov))


    def calculate_basic(self):
        z_height=int(self.z_list.shape[0]/self.N)
        x_height=int(self.Cx_list.shape[0]/self.N)

        for i in range(self.N):
            z=self.z_list[(z_height*i):(z_height*i+z_height),:]
            z_t=self._z_t_list[:,(z_height*i):(z_height*i+z_height)]
            x=self.Cx_list[(x_height*i):(x_height*i+x_height),:]
            y = self.Cy_list[(x_height * i):(x_height * i + x_height), :]
            if i==0:
                temp_xz= (z @ x).transpose()
                temp_zy = (z @ y).transpose()
            else:
                temp_xz += (z @ x).transpose()
                temp_zy += (z @ y).transpose()
        return (temp_xz, temp_zy)

    def calculate_W(self, H):
        #W1 = (1.0 / N) * sum_product2([z_list, H1, _z_t_list], [(N, 1), (1, 1), (1, N)])
        z_height = int(self.z_list.shape[0] / self.N)

        for i in range(self.N):
            z = self.z_list[(z_height * i):(z_height * i + z_height), :]
            z_t = self._z_t_list[:, (z_height * i):(z_height * i + z_height)]

            if i == 0:
                temp_W = z @ H @ z_t

            else:
                temp_W += z @ H @ z_t

        return temp_W


    def calculate_residual(self, y_list, x_list, beta):

        tbr=y_list - x_list @ beta

        return (tbr)

    def vcov(self, step: int):
        # report robust vcov only
        step_1 = self.result_list[0]
        step_2 = self.result_list[1]

        if step >= 2:
            the_step=self.result_list[step-1]
            previous_step=self.result_list[step-2]
            M2=the_step.M
            _M2_XZ_W2=the_step._M_XZ_W
            _W2_inv=the_step.W_inv
            zs2=the_step.zs
            vcov_step1=previous_step.vcov
            residual1=previous_step.residual
            return Windmeijer(M2, _M2_XZ_W2, _W2_inv, zs2,
                                               vcov_step1, self.Cx_list, self.z_list, residual1, self.N)
        elif step==1:
            _M_XZ_W=step_1._M_XZ_W
            W2=step_2.W
            return self.N * (_M_XZ_W @ W2 @ _M_XZ_W.transpose())

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
        width = z_list.shape[1]

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
        xy_height=int(Cy_list.shape[0]/N)
        z_height=int(z_list.shape[0]/N)

        na_list = []
        num_NA=0
        #num_NA_diff = 0
        #num_NA_level=0


        # max_observations_per_group = int((np.size(Cy_list[0]) + 1) / 2)

        for i in range(N):
            x = Cx_list[(i*xy_height):(i*xy_height+xy_height),:]
            y = Cy_list[(i*xy_height):(i*xy_height+xy_height),:]
            z = z_list[(i*z_height):(i*z_height+z_height),:]
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


    def generate_summary(self, step):
        step1 = self.result_list[0]
        step2 = self.result_list[1]
        if step==1 or step==2:
            _W2_inv=step2.W_inv
            zs=step2.zs

            self.hansen = tests.hansen_overid(_W2_inv, self.N, zs, self.num_instru, \
                                          self.Cx_list.shape[1])
            
            
        else:
            current_step=self.result_list[step-1]
            _W2_inv = current_step.W_inv
            zs = current_step.zs
            self.hansen = tests.hansen_overid(_W2_inv, self.N, zs, self.num_instru, \
                                              self.Cx_list.shape[1])

        self.AR_list = tests.AR_test(self, step, 2)

        if self.steps==2:
            str_steps='two-step '
        elif self.steps==1:
            str_steps='one-step '
        else:
            str_steps=str(self.steps) + '-step '
        if self.level:
            str_gmm='system GMM'
        else:
            str_gmm='difference GMM'

        print('Dynamic panel-data estimation, ' + str_steps + str_gmm)
        print(self.basic_information())
        print(self.regression_table(step))
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

    def regression_table(self, step):
        regression_result=self.result_list[step-1]
        beta = regression_result.beta
        std_err = regression_result.std_err

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