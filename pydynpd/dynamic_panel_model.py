import numpy as np

from pydynpd.info import df_info, options_info
from pydynpd.instruments import instruments
from pydynpd.panel_data import panel_data
from pydynpd.variable import regular_variable


class data_table(object):
    def __init__(self, dat: np.ndarray, unit_height: int):
        self.dat = dat
        self.width = dat.shape[1]
        self.height = dat.shape[0]
        self.unit_height = unit_height


class dynamic_panel_model(object):
    def __init__(self, pdata: panel_data, variables: dict, options: options_info):
        self.pdata = pdata
        self.T = self.pdata.T
        self.N = self.pdata.N
        self.variables = variables.copy()
        self.options = options
        method = 'fd'

        max_lag, first_diff_index, first_level_index, last_index = self.get_info(variables, 'fd', self.T)

        self.df_information = df_info(N=self.N, T=self.T, ids=self.pdata.ids, max_lag=max_lag,
                                      first_diff_index=first_diff_index, first_level_index=first_level_index,
                                      last_index=last_index)

        self.step_results = []

        if first_diff_index + 2 > last_index:  # to do: change 3 to something rated to AR(p)
            raise Exception("Not enough periods to run the model")
        else:
            if options.timedumm:
                self.update_time_dummies(first_diff_index, last_index)
            self.prepare_data()

    def update_time_dummies(self, first_diff_index, last_index):

        # for t in range(first_index-1, last_index+1):
        #     var_name=self.pdata.col_timedumm[t]
        #     if t==first_index-1:
        #         if self.options.level:
        #             new_var = regular_variable(var_name, 0)
        #             self.variables['iv'].append(new_var)
        #     else:
        #         new_var = regular_variable(var_name, 0)
        #         self.variables['dep_indep'].append(new_var)
        #         self.variables['iv'].append(new_var)

        for var_name in self.pdata.col_timedumm[first_diff_index:(last_index + 1)]:
            new_var = regular_variable(var_name, 0)
            self.variables['dep_indep'].append(new_var)
            self.variables['iv'].append(new_var)

    def prepare_data(self):
        gmm_tables = self.get_gmm_table_dict(self.variables, self.options.level)
        xy_tables = self.get_xy_table_dict(self.variables)
        self.final_xy_tables = self.get_final_xy_tables(xy_tables, self.options.level)
        z = instruments(self.variables, gmm_tables, self.df_information, self.options)
        self.z_information = z.z_information
        self.z_list = z.z_table
        self.num_obs = self.prepare_reg()
        self._z_t_list = self.z_list.transpose()

    def get_info(self, variables, method, T):
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

        if method == 'fd':
            last_index = T - 1  # zero based
            first_diff_index = max(max_lag + 1, max_Dgmm_minlag)
            first_level_index = max(max_lag, max_Lgmm_minlag)

        else:
            last_index = T - 2
            first_diff_index = max_lag - 1

        return (max_lag, first_diff_index, first_level_index, last_index)

    def get_gmm_table_dict(self, variables, level):
        gmm_dict = {}

        Dgmm_table = self.gen_table(self.pdata.data, variables['Dgmm'], [])
        iv_table = self.gen_table(self.pdata.data, variables['iv'], [])
        Delta_iv_table = self.gen_table(self.pdata.fd_data, variables['iv'], [])

        gmm_dict['Dgmm'] = Dgmm_table
        gmm_dict['iv'] = iv_table
        gmm_dict['Div'] = Delta_iv_table

        if level:  # sys-GMM
            Lgmm_table = self.gen_table(self.pdata.fd_data, variables['Lgmm'], [])
            gmm_dict['Lgmm'] = Lgmm_table

        return (gmm_dict)

    def get_xy_table_dict(self, variables: dict):
        xy_tables = {}

        num_var = len(variables['dep_indep'])
        cut = [self.df_information.first_level_index, self.df_information.last_index]
        Dcut = [self.df_information.first_diff_index, self.df_information.last_index]

        y_table = self.gen_table(self.pdata.data, variables['dep_indep'][0:1], cut)
        x_table = self.gen_table(self.pdata.data, variables['dep_indep'][1:num_var], cut)

        Dy_table = self.gen_table(self.pdata.fd_data, variables['dep_indep'][0:1], Dcut)
        Dx_table = self.gen_table(self.pdata.fd_data, variables['dep_indep'][1:num_var], Dcut)

        xy_tables['x'] = x_table
        xy_tables['y'] = y_table
        xy_tables['Dx'] = Dx_table
        xy_tables['Dy'] = Dy_table

        return xy_tables

    def get_final_xy_tables(self, xy_tables, level):
        final_xy_tables = {}

        N = self.pdata.N

        Dx_dat = xy_tables['Dx'].dat
        Dy_dat = xy_tables['Dy'].dat
        x_dat = xy_tables['x'].dat
        y_dat = xy_tables['y'].dat

        Dx_height = xy_tables['Dx'].unit_height  # int(Dx_dat.shape[0] / N)
        x_height = xy_tables['x'].unit_height  # int(x_dat.shape[0] / N)

        height_total = Dx_height + x_height

        if level:  # sys-GMM
            width = x_dat.shape[1]
            Cx_dat = np.empty((height_total * N, width + 1), dtype=np.float64)
            Cy_dat = np.empty((height_total * N, 1), dtype=np.float64)
            # Cx=np.empty()
            for i in range(N):  # , nogil=True):  #df_information.N
                Dy_i = Dy_dat[(i * Dx_height):(i * Dx_height + Dx_height), :]
                y_i = y_dat[(i * x_height):(i * x_height + x_height), :]

                temp_y = Cy_dat[(height_total * i):(height_total * (i + 1)), :]
                temp_y[0:Dx_height, 0] = Dy_i[0:Dx_height, 0]
                temp_y[Dx_height:height_total, 0] = y_i[0:x_height, 0]

                Dx_i = Dx_dat[(i * Dx_height):(i * Dx_height + Dx_height), :]
                x_i = x_dat[(i * x_height):(i * x_height + x_height), :]

                temp_x = Cx_dat[(height_total * i):(height_total * (i + 1)), :]
                temp_x[0:Dx_height, 0:width] = Dx_i[0:Dx_height, 0:width]
                temp_x[Dx_height:height_total, 0:width] = x_i[0:x_height, 0:width]
                temp_x[0:Dx_height, width] = 0
                temp_x[Dx_height:height_total, width] = 1

            Cy = data_table(Cy_dat, height_total)
            Cx = data_table(Cx_dat, height_total)

        else:  # diff-GMM
            Cx = xy_tables['Dx']
            Cy = xy_tables['Dy']

        final_xy_tables['Cy'] = Cy
        final_xy_tables['Cx'] = Cx

        return (final_xy_tables)

    def gen_table(self, ori_data: np.ndarray, variable_list, cut: list):
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

        if len(cut) == 0:
            start_row = 0
            end_row = T - 1
        else:
            start_row = cut[0]
            end_row = cut[1]

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
                    if len(cut) == 0:
                        tbr_i[0:var.lag, col] = np.NaN
                        tbr_i[var.lag:(end_row + 1), col] = ori_i[0:(end_row + 1 - var.lag), which_col[j]]
                    else:
                        tbr_i[:, col] = ori_i[(start_row - var.lag):(end_row + 1 - var.lag), which_col[j]]

                col += 1

        new_table = data_table(tbr, height)
        return (new_table)

    def prepare_reg(self):

        z_list = self.z_list
        # num_instru = model.z_information.num_instr
        Cx = self.final_xy_tables['Cx']
        Cy = self.final_xy_tables['Cy']
        Cx_dat = Cx.dat
        Cy_dat = Cy.dat

        N = self.N
        xy_height = Cy.unit_height
        z_height = int(z_list.shape[0] / N)

        na_list = []
        num_NA = 0

        for i in range(N):
            x = Cx_dat[(i * xy_height):(i * xy_height + xy_height), :]
            y = Cy_dat[(i * xy_height):(i * xy_height + xy_height), :]
            z = z_list[(i * z_height):(i * z_height + z_height), :]
            row_if_nan = np.logical_or(np.isnan(x).any(axis=1), np.isnan(y).any(axis=1))
            # num_NA+=np.count_nonzero(row_if_nan[range(0,int((np.size(y)+1)/2))])
            # num_NA2 += np.count_nonzero(row_if_nan[range(max_observations_per_group,??? )])
            num_NA += np.count_nonzero(row_if_nan)
            if self.options.level:
                num_NA -= np.count_nonzero(row_if_nan[range(0, self.z_information.diff_width)])

            na_list.append(row_if_nan)
            for j in range(0, len(row_if_nan)):
                if row_if_nan[j] == True:
                    x[j, :] = 0
                    y[j, :] = 0
                    z[:, j] = 0
                    # n_obs = n_obs - 1
        # return(Cx_dat, Cy_dat, z_list)

        if self.options.level:
            return (self.z_information.level_width * N - num_NA)
        else:
            return (self.z_information.diff_width * N - num_NA)
