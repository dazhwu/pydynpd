import numpy as np

from pydynpd.info import df_info, options_info
from pydynpd.panel_data import panel_data
from pydynpd.instruments import instruments
from pydynpd.variable import regular_variable


class dynamic_panel_model(object):
    def __init__(self, pdata: panel_data, variables: dict, options: options_info):
        self.pdata = pdata
        self.T = self.pdata.T
        self.N = self.pdata.N
        self.variables = variables.copy()
        self.options = options
        method='fd'

        max_lag, first_index, last_index = self.get_info(variables, 'fd', self.T)

        self.df_information = df_info(N=self.N, T=self.T, ids=self.pdata.ids, max_lag=max_lag,
                                      first_index=first_index, last_index=last_index)

        
        self.step_results=[]
        self.results={}

        if first_index + 2 > last_index:  # to do: change 3 to something rated to AR(p)
            raise Exception("Not enough periods to run the model")
        else:
            if options.timedumm:
                self.update_time_dummies(first_index, last_index)
            self.prepare_data()
    
    def update_time_dummies(self, first_index, last_index):       
        
        for var_name in self.pdata.col_timedumm[first_index:(last_index+1)]:
            new_var = regular_variable(var_name, 0)
            self.variables['dep_indep'].append(new_var)
            self.variables['iv'].append(new_var)


    def prepare_data(self):
        gmm_tables = self.get_gmm_table_dict(self.variables, self.options.level)
        self.xy_tables = self.get_xy_table_dict(self.variables)
        self.final_xy_tables = self.get_final_xy_tables(self.options.level)
        z = instruments(self.variables, gmm_tables,self.df_information,self.options )
        self.z_information=z.z_information
        self.z_list =z.z_list
        self.num_obs = self.prepare_reg()
        self._z_t_list=self.z_list.transpose()

    def get_info(self,variables, method, T):
        max_lag = 0
        max_gmm_minlag = 0
        for var in variables['dep_indep'] + variables['iv']:
            if var.lag > max_lag:
                max_lag = var.lag

        for var in variables['gmm']:
            if var.min_lag > max_gmm_minlag:
                max_gmm_minlag = var.min_lag

        if method == 'fd':
            last_index = T - 1  # zero based
            first_index = max(max_lag + 1, max_gmm_minlag)
        else:
            last_index = T - 2
            first_index = max_lag - 1

        return (max_lag, first_index, last_index)

    def get_gmm_table_dict(self, variables, level):
        iv_table = self.gen_table(self.pdata.data, variables['iv'], [])
        gmm_table = self.gen_table(self.pdata.data, variables['gmm'], [])
        Div_table = self.gen_table(self.pdata.fd_data, variables['iv'], [])

        gmm_dict = {}

        if level:  # sys-GMM
            Dgmm_table = self.gen_table(self.pdata.fd_data, variables['gmm'], [])
            gmm_dict['Dgmm'] = Dgmm_table

        gmm_dict['gmm'] = gmm_table
        gmm_dict['iv'] = iv_table

        gmm_dict['Div'] = Div_table

        return (gmm_dict)

    def get_xy_table_dict(self, variables: dict):
        xy_tables = {}

        num_var = len(variables['dep_indep'])
        cut = [self.df_information.first_index - 1, self.df_information.last_index]
        Dcut = [self.df_information.first_index, self.df_information.last_index]

        y_table = self.gen_table(self.pdata.data, variables['dep_indep'][0:1], cut)
        x_table = self.gen_table(self.pdata.data, variables['dep_indep'][1:num_var], cut)

        Dy_table = self.gen_table(self.pdata.fd_data, variables['dep_indep'][0:1], Dcut)
        Dx_table = self.gen_table(self.pdata.fd_data, variables['dep_indep'][1:num_var], Dcut)

        xy_tables['x'] = x_table
        xy_tables['y'] = y_table
        xy_tables['Dx'] = Dx_table
        xy_tables['Dy'] = Dy_table

        return xy_tables

    def get_final_xy_tables(self, level):
        final_xy_tables = {}

        N = self.pdata.N

        Dx_list = self.xy_tables['Dx']
        Dy_list = self.xy_tables['Dy']
        x_list = self.xy_tables['x']
        y_list = self.xy_tables['y']

        Dx_height = int(Dx_list.shape[0] / N)
        x_height = int(x_list.shape[0] / N)

        height_total = Dx_height + x_height

        if level:  # sys-GMM
            width = x_list.shape[1]
            Cx_list = np.empty((height_total * N, width + 1), dtype=np.float64)
            Cy_list = np.empty((height_total * N, 1), dtype=np.float64)
            # Cx=np.empty()
            for i in range(N):  # , nogil=True):  #df_information.N
                Dy_i = Dy_list[(i * Dx_height):(i * Dx_height + Dx_height), :]
                y_i = y_list[(i * x_height):(i * x_height + x_height), :]

                temp_y = Cy_list[(height_total * i):(height_total * (i + 1)), :]
                temp_y[0:Dx_height, 0] = Dy_i[0:Dx_height, 0]
                temp_y[Dx_height:height_total, 0] = y_i[0:x_height, 0]

                Dx_i = Dx_list[(i * Dx_height):(i * Dx_height + Dx_height), :]
                x_i = x_list[(i * x_height):(i * x_height + x_height), :]

                temp_x = Cx_list[(height_total * i):(height_total * (i + 1)), :]
                temp_x[0:Dx_height, 0:width] = Dx_i[0:Dx_height, 0:width]
                temp_x[Dx_height:height_total, 0:width] = x_i[0:x_height, 0:width]
                temp_x[0:Dx_height, width] = 0
                temp_x[Dx_height:height_total, width] = 1

        else:  # diff-GMM
            Cx_list = Dx_list
            Cy_list = Dy_list

        final_xy_tables['Cy'] = Cy_list
        final_xy_tables['Cx'] = Cx_list

        return (final_xy_tables)

    def gen_table(self, ori_data: np.ndarray, variable_list, cut: list):
        num_variables = len(variable_list)
        list_cols = self.pdata.ids.copy()
        N = self.pdata.N
        T = self.pdata.T

        for var in variable_list:
            if var.name not in list_cols:
                list_cols.append(var.name)

        
        df_cols=self.pdata.cols + self.pdata.col_timedumm
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
            ori_i = ori_data[(i * T ):(i * T + T), :]
            tbr_i = tbr[(i * height):(i * height + height), :]
            for j in range(num_variables):
                var = variable_list[j]
                if var.lag == 0:
                    tbr_i[:, col] = ori_i[start_row:(end_row+1), which_col[j]]
                else: 
                    if len(cut)==0:
                        tbr_i[0:var.log,col]=np.NaN
                        tbr_i[var.lag:(end_row+1),col]=ori_i[0:(end_row+1-var.lag),which_col[j]]
                    else:
                        tbr_i[:,col]=ori_i[(start_row-var.lag):(end_row+1-var.lag),which_col[j]]   
                    
                col += 1

        return (tbr)
    
    def prepare_reg(self):

        z_list=self.z_list
        #num_instru = model.z_information.num_instr
        Cx_list = self.final_xy_tables['Cx']
        Cy_list = self.final_xy_tables['Cy']

        N = self.N
        xy_height = int(Cy_list.shape[0] / N)
        z_height = int(z_list.shape[0] / N)

        na_list = []
        num_NA = 0

        for i in range(N):
            x = Cx_list[(i * xy_height):(i * xy_height + xy_height), :]
            y = Cy_list[(i * xy_height):(i * xy_height + xy_height), :]
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
        # return(Cx_list, Cy_list, z_list)

        if self.options.level:
            return (self.z_information.level_width * N - num_NA)
        else:
            return (self.z_information.diff_width * N - num_NA)

    def form_results(self):
        step=len(self.step_results)
        the_list=self.step_results[step-1]
        self.results['beta']=the_list.beta
        self.results['std_err']=the_list.std_err
        self.results['vcov']=the_list.vcov
        self.results['W']=the_list.W