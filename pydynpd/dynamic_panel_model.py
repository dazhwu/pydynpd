import math

import numpy as np
import pandas as pd
import scipy

from pydynpd.common_functions import get_first_diff_table, get_fod_table
from pydynpd.info import df_info, options_info
from pydynpd.instruments import instruments
from pydynpd.panel_data import panel_data
from pydynpd.variable import regular_variable


class dynamic_panel_model(object):
    def __init__(self, pdata: panel_data, variables: dict, options: options_info, command_str: str, part_2: str,
                 part_3: str):
        self.name = ''
        self.pdata = pdata
        self.T = self.pdata.T
        self.N = self.pdata.N
        self.variables = variables.copy()
        self.options = options

        self.command_str = command_str + '|' + part_2
        if part_3 != '':
            self.command_str += '|' + part_3

        self.get_info(variables, self.N, self.T)

        self.step_results = []

        if options.timedumm:
            self.update_time_dummies(self.df_information.first_diff_index, self.df_information.last_diff_index)

        self.prepare_data()

    def update_time_dummies(self, first_diff_index, last_index):

        for var_name in self.pdata.col_timedumm[first_diff_index:(last_index + 1)]:
            new_var = regular_variable(var_name, 0)
            self.variables['dep_indep'].append(new_var)
            self.variables['iv'].append(new_var)

    def prepare_data(self):
        gmm_tables = self.get_gmm_table_dict(self.variables, self.options.level)
        xy_tables = self.get_xy_table_dict(self.variables)
        self.final_xy_tables = self.get_final_xy_tables(xy_tables, self.options.level, self.options.transformation)
        z = instruments(self.variables, gmm_tables, self.df_information, self.options)
        self.z_information = z.z_information
        self.z_list = z.z_table
        self.num_obs, self.max_obs, self.min_obs, self.avg_obs = self.prepare_reg()
        self._z_t_list = self.z_list.transpose()

    def get_info(self, variables, N, T):
        max_lag = 0
        max_Dgmm_minlag = 0
        max_Lgmm_minlag = 0
        for var in variables['dep_indep'] + variables['iv']:
            if var.lag > max_lag:
                max_lag = var.lag

        for var in variables['Dgmm']:
            if var.min_lag > max_Dgmm_minlag:
                max_Dgmm_minlag = var.min_lag

        for var in variables['Lgmm']:
            if var.min_lag > max_Lgmm_minlag:
                max_Lgmm_minlag = var.min_lag

        last_level_index = T - 1
        last_diff_index = T - 1

        first_level_index = max(max_lag, max_Lgmm_minlag)

        if self.options.level and self.options.transformation == 'fod':
            first_diff_index = first_level_index
        else:
            first_diff_index = first_level_index + 1  # max(max_lag + 1, max_Dgmm_minlag)

        if first_diff_index + 2 > last_diff_index:  # to do: change 3 to something rated to AR(p)
            raise Exception("Not enough periods to run the model")

        self.df_information = df_info(N=self.N, T=self.T, ids=self.pdata.ids, max_lag=max_lag,
                                      first_diff_index=first_diff_index, last_diff_index=last_diff_index,
                                      first_level_index=first_level_index,
                                      last_level_index=last_level_index)

    def get_gmm_table_dict(self, variables, level):
        gmm_dict = {}

        Dgmm_table = self.gen_table(self.pdata.data, variables['Dgmm'])
        iv_table = self.gen_table(self.pdata.data, variables['iv'])
        if (self.options.transformation=='fod') & (level==False):
            Delta_iv_table=self.gen_table(self.pdata.fod_data, variables['iv'])
        else:
            Delta_iv_table = self.gen_table(self.pdata.fd_data, variables['iv'])

        gmm_dict['Dgmm'] = Dgmm_table
        gmm_dict['iv'] = iv_table

        gmm_dict['Div'] = Delta_iv_table

        if level:  # sys-GMM
            Lgmm_table = self.gen_table(self.pdata.fd_data, variables['Lgmm'])
            gmm_dict['Lgmm'] = Lgmm_table

        return (gmm_dict)

    def get_xy_table_dict(self, variables: dict):
        xy_tables = {}

        num_var = len(variables['dep_indep'])
        ori_y_table = self.gen_table(self.pdata.data, variables['dep_indep'][0:1])
        ori_x_table = self.gen_table(self.pdata.data, variables['dep_indep'][1:num_var])

        xy_tables['x'] = ori_x_table
        xy_tables['y'] = ori_y_table

        ori_Diff_y_table = get_first_diff_table(ori_y_table, self.N)
        ori_Diff_x_table = get_first_diff_table(ori_x_table, self.N)

        if self.options.transformation == 'fd':
            xy_tables['Dy'] = ori_Diff_y_table
            xy_tables['Dx'] = ori_Diff_x_table

        else:
            ori_Fod_y_table = get_fod_table(ori_y_table, self.N)
            ori_Fod_x_table = get_fod_table(ori_x_table, self.N)

            xy_tables['Dy'] = ori_Fod_y_table
            xy_tables['Dx'] = ori_Fod_x_table

            xy_tables['Diff_y'] = ori_Diff_y_table
            xy_tables['Diff_x'] = ori_Diff_x_table

        return xy_tables

    def get_final_xy_tables(self, xy_tables, level, transformation):
        final_xy_tables = {}
        N = self.pdata.N

        Dcut = [self.df_information.first_diff_index, self.df_information.last_diff_index]
        cut = [self.df_information.first_level_index, self.df_information.last_level_index]

        Dcut_height = Dcut[1] - Dcut[0] + 1
        cut_height = cut[1] - cut[0] + 1

        if level:  # sys-GMM
            Cy, Cx = self.get_final_xy_systemGMM(xy_tables, transformation, cut, Dcut, cut_height, Dcut_height)
        else:  # diff-GMM
            Cy, Cx = self.get_final_xy_diffGMM(xy_tables['Dy'], xy_tables['Dx'], Dcut, Dcut_height)

        final_xy_tables['Cy'] = Cy
        final_xy_tables['Cx'] = Cx

        if self.options.transformation == 'fod':
            Diff_y, Diff_x = self.get_final_xy_diffGMM(xy_tables['Diff_y'], xy_tables['Diff_x'], Dcut, Dcut_height)

            if self.options.level:
                height = Diff_y.shape[0]
                zeros = np.zeros((height, 1), dtype=np.float64)
                # Diff_y=np.hstack((Diff_y, zeros))
                # zero_xs = zeros.copy()
                Diff_x = np.hstack((Diff_x, zeros))

            final_xy_tables['Diff_y'] = Diff_y
            final_xy_tables['Diff_x'] = Diff_x

        return (final_xy_tables)

    def get_final_xy_systemGMM(self, xy_tables, transformation, cut, Dcut, cut_height, Dcut_height):
        Dx = xy_tables['Dx']
        Dy = xy_tables['Dy']
        x = xy_tables['x']
        y = xy_tables['y']

        N = self.N

        height_total = Dcut_height + cut_height

        width = x.shape[1]
        Dx_height = int(Dx.shape[0] / N)
        x_height = int(x.shape[0] / N)

        Cx = np.empty((height_total * N, width + 1), dtype=np.float64)
        Cy = np.empty((height_total * N, 1), dtype=np.float64)
        # Cx=np.empty()
        for i in range(N):  # , nogil=True):  #df_information.N
            Dy_i = Dy[(i * Dx_height):(i * Dx_height + Dx_height), :]
            y_i = y[(i * x_height):(i * x_height + x_height), :]

            temp_y = Cy[(height_total * i):(height_total * (i + 1)), :]
            temp_y[0:Dcut_height, 0] = Dy_i[Dcut[0]:(Dcut[1] + 1), 0]
            temp_y[Dcut_height:height_total, 0] = y_i[cut[0]:(cut[1] + 1), 0]

            if transformation == 'fod':
                temp_y[0:1, 0] = np.NaN

            Dx_i = Dx[(i * Dx_height):(i * Dx_height + Dx_height), :]
            x_i = x[(i * x_height):(i * x_height + x_height), :]

            temp_x = Cx[(height_total * i):(height_total * (i + 1)), :]
            temp_x[0:Dcut_height, 0:width] = Dx_i[Dcut[0]:(Dcut[1] + 1), 0:width]
            temp_x[Dcut_height:height_total, 0:width] = x_i[cut[0]:(cut[1] + 1), 0:width]
            temp_x[0:Dcut_height, width] = 0
            temp_x[Dcut_height:height_total, width] = 1

        return (Cy, Cx)

    def get_final_xy_diffGMM(self, Dy, Dx, Dcut, Dcut_height):

        N = self.N

        Dx_height = int(Dx.shape[0] / N)

        width = Dx.shape[1]

        Cx = np.empty((Dcut_height * N, width), dtype=np.float64)
        Cy = np.empty((Dcut_height * N, 1), dtype=np.float64)

        for i in range(N):  # , nogil=True):  #df_information.N
            Dy_i = Dy[(i * Dx_height):(i * Dx_height + Dx_height), :]

            temp_y = Cy[(Dcut_height * i):(Dcut_height * (i + 1)), :]
            temp_y[0:Dcut_height, 0] = Dy_i[Dcut[0]:(Dcut[1] + 1), 0]

            Dx_i = Dx[(i * Dx_height):(i * Dx_height + Dx_height), :]

            temp_x = Cx[(Dcut_height * i):(Dcut_height * (i + 1)), :]
            temp_x[0:Dcut_height, 0:width] = Dx_i[Dcut[0]:(Dcut[1] + 1), 0:width]

        return (Cy, Cx)

    def gen_table(self, ori_data: np.ndarray, variable_list):
        num_variables = len(variable_list)
        list_cols = self.pdata.ids.copy()
        N = self.pdata.N
        T = self.pdata.T

        for var in variable_list:
            if var.name not in list_cols:
                list_cols.append(var.name)

        df_cols = self.pdata.cols + self.pdata.col_timedumm
        variable_names = [var.name for var in variable_list]
        which_col = [df_cols.index(var_name) for var_name in variable_names]

        start_row = 0
        end_row = T - 1

        height = end_row - start_row + 1
        tbr = np.empty((height * N, num_variables),
                       dtype='float64')  # ori_arr[(info.first_index - 1):(info.last_index + 1), :]

        for i in range(N):
            col = 0
            ori_i = ori_data[(i * T):(i * T + T), :]
            tbr_i = tbr[(i * height):(i * height + height), :]

            for j in range(num_variables):

                var = variable_list[j]
                if var.lag == 0:
                    tbr_i[:, col] = ori_i[start_row:(end_row + 1), which_col[j]]

                else:
                    tbr_i[0:var.lag, col] = np.NaN
                    tbr_i[var.lag:(end_row + 1), col] = ori_i[0:(end_row + 1 - var.lag), which_col[j]]

                col += 1

        return tbr

    def calculate_MMSC_LU(self):
        self.MMSC_LU = {}
        log_n = math.log(self.num_obs)
        dif = self.z_information.num_instr - (len(self.variables['dep_indep']) - 1)
        self.MMSC_LU["bic"] = self.hansen.test_value - (dif) * log_n
        self.MMSC_LU["hqic"] = self.hansen.test_value - dif * math.log(log_n) * 2.1
        self.MMSC_LU["aic"] = self.hansen.test_value - (dif) * 2

    def form_regression_table(self):

        num_indeps = len(self.variables['dep_indep']) - 1

        var_names = []
        for i in range(1, num_indeps + 1):
            var_name = self.variables['dep_indep'][i].name
            var_lag = self.variables['dep_indep'][i].lag
            if (var_lag) >= 1:
                var_name = 'L' + str(var_lag) + '.' + var_name
            var_names.append(var_name)

        if self.options.level:
            var_names.append('_con')
            num_indeps += 1

        num_steps = len(self.step_results)
        if self.options.steps <= 2:
            the_result = self.step_results[self.options.steps - 1]
        else:
            the_result = self.step_results[num_steps - 1]

        coeff = the_result.beta[:, 0]
        std_err = the_result.std_err
        z_value = [coeff[i] / std_err[i] for i in range(num_indeps)]
        p_value = [scipy.stats.norm.sf(abs(z)) * 2 for z in z_value]
        sig = ['***' if p <= 0.001 else ('**' if p <= 0.01 else ('*' if p <= 0.05 else ' ')) for p in p_value]

        self.regression_table = pd.DataFrame(list(zip(var_names, coeff, std_err, z_value, p_value, sig)),
                                             columns=['variable', 'coefficient', 'std_err', 'z_value', 'p_value',
                                                      'sig'])

    def prepare_reg_fod(self):

        Diff_x = self.final_xy_tables['Diff_x']
        Diff_y = self.final_xy_tables['Diff_y']
        N = self.N
        xy_height = int(Diff_x.shape[0] / N)
        row_if_nan = np.logical_or(np.isnan(Diff_y).any(axis=1), np.isnan(Diff_x).any(axis=1))

        for i in range(len(row_if_nan)):
            if row_if_nan[i]:
                Diff_x[i, :] = 0
                Diff_y[i, :] = 0

    def prepare_reg(self):

        z_list = self.z_list

        Cx = self.final_xy_tables['Cx']
        Cy = self.final_xy_tables['Cy']

        N = self.N
        xy_height = int(Cy.shape[0] / N)
        z_height = int(z_list.shape[0] / N)

        na_list = []
        num_NA = 0
        total = 0
        max = 0
        min = 0
        if self.options.transformation == 'fod':
            self.prepare_reg_fod()
        for i in range(N):
            x = Cx[(i * xy_height):(i * xy_height + xy_height), :]
            y = Cy[(i * xy_height):(i * xy_height + xy_height), :]
            z = z_list[(i * z_height):(i * z_height + z_height), :]
            row_if_nan = np.logical_or(np.isnan(x).any(axis=1), np.isnan(y).any(axis=1))

            if self.options.level:
                temp = np.count_nonzero(row_if_nan[range(self.z_information.diff_width, self.z_information.width)])
            else:
                temp = np.count_nonzero(row_if_nan)
            num_NA += temp

            if temp > max:
                max = temp
            if temp < min:
                min = temp

            na_list.append(row_if_nan)
            for j in range(0, len(row_if_nan)):
                if row_if_nan[j] == True:
                    x[j, :] = 0
                    y[j, :] = 0
                    z[:, j] = 0
                    # if self.options.transformation=='fod':
                    #
                    #     Diff_x[j+i*self.z_information.diff_width,:]=0
                    #     Diff_y[j + i * self.z_information.diff_width, :] = 0

        if self.options.level:
            width = self.z_information.level_width
        else:
            width = self.z_information.diff_width

        nobs = width * N - num_NA
        max_obs = width - min
        min_obs = width - max
        avg_obs = width - (num_NA * 1.0) / N

        return (nobs, max_obs, min_obs, avg_obs)
